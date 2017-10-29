package org.orangepalantir.genericpca;

import lightgraph.DataSet;
import lightgraph.Graph;
import org.orangepalantir.genericpca.CoefficientKmeansND;
import org.orangepalantir.genericpca.IndexedCoefficient;

import javax.swing.JTextField;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Created by msmith on 19.10.17.
 */
public class ScanKmeans {
    static final int MAXRESULTS = 50;
    PriorityQueue<Result> results = new PriorityQueue<>(MAXRESULTS + 1);
    List<Condition> conditions = new ArrayList<>();
    int counted = 0;
    int finished = 0;
    List<List<IndexedCoefficient>> coefficients;
    static final int THREADS = 4;
    ExecutorService service;
    public ScanKmeans(List<List<IndexedCoefficient>> coefficients){
        this.coefficients = coefficients;
    }

    private void generateConditions(){
        int total= 512;
        int count = total;
        int alt = total - 1;
        int small = alt - 1;
        int top = coefficients.get(0).size() - 1;

        int n = 1;




        while(count>=2){
            int a = top - count;
            int b = top - alt;
            int c = top - small;
            int[] indexes = {a, b, c};
            small--;
            if(small<0){
                alt--;
                if(alt<1){
                    count--;
                    if(count<2){
                        break;
                    }
                    alt = count -1;
                }

                small = alt-1;
            }


            for(int k = 0; k<n; k++){
                int km = 6 + k;
                conditions.add(new Condition(indexes, km));
            }
        }
        System.out.println("generated");

        Collections.shuffle(conditions);

    }

    synchronized public void postResult(Worker w, Result r){
        results.add(r);
        while(results.size()>MAXRESULTS){
            results.poll();
        }

        finished++;
        if(finished==conditions.size()){
            System.out.println("finished");
            service.shutdownNow();
            notifyAll();
            return;
        }


        if(counted<conditions.size()){
            w.post(conditions.get(counted));
            counted++;
        }
        if(counted%10000==0){
            logHighest();
        }
    }
    public void logHighest(){
        try(BufferedWriter writer = Files.newBufferedWriter(
                Paths.get("bestKs.txt"),
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING
        )) {
            for(Result result: results){
                int[] dexes = result.indexes;
                writer.write(String.format("%d\t%d\t%d\t%d\t%f\n", dexes[0], dexes[1], dexes[2], result.ks, result.variance));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    synchronized public void start(){
        System.out.println("generating");
        generateConditions();
        System.out.println("conditions generated");
        service = Executors.newFixedThreadPool(THREADS);
        for(int i = 0; i<THREADS; i++){
            Worker worker = new Worker(coefficients);
            service.submit(worker);
            worker.post(conditions.get(counted));
            counted++;
        }

    }
    public void join(){
        while(finished<conditions.size()){
            synchronized (this){
                try {
                    wait();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }
        logHighest();
    }


    public static void main(String[] args) throws IOException {
        List<List<IndexedCoefficient>> coefficients = IndexedCoefficient.readCoefficients(Paths.get(args[0]));
        ScanKmeans skm = new ScanKmeans(coefficients);

        skm.start();
        skm.join();
        //graph.show(true, "Variance for different Km values");
    }
    static class Condition{
        int[] indexes;
        int k;
        public Condition(int[] indexes, int k){
            this.indexes = indexes;
            this.k = k;
        }
    }

    class Worker implements Runnable{
        BlockingQueue<Condition> queue = new LinkedBlockingQueue<>(1);
        CoefficientKmeansND kmeans;
        int attempts = 5;
        public Worker(List<List<IndexedCoefficient>> coefficients){
            kmeans = new CoefficientKmeansND();
            kmeans.setInput(coefficients);

        }

        public void post(Condition condition){
            queue.add(condition);
        }

        @Override
        public void run() {
            System.out.println("starting");
            while(!Thread.interrupted()){
                try {
                    Condition condition = queue.take();
                    kmeans.ks = condition.k;
                    double s = Double.MAX_VALUE;
                    for(int i = 0; i<attempts; i++){
                        double s2 = kmeans.calculate(condition.indexes);
                        s = s2<s?s2:s;
                    }
                    Result result = new Result(condition.indexes, s, condition.k);
                    postResult(this, result);
                } catch (InterruptedException e) {
                    return;
                }
            }
        }

    }

    static class Result implements Comparable<Result>{
        int[] indexes;
        double variance;
        int ks;
        public Result(int[] is, double v, int k){
            indexes = is;
            variance = v;
            ks = k;
        }

        @Override
        public int compareTo(Result o) {
            return -Double.compare(variance, o.variance);
        }
        @Override
        public String toString(){
            return variance + " : " + Arrays.toString(indexes) + " " + ks;
        }
    }

}
