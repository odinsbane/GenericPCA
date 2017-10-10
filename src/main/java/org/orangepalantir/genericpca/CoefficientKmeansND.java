package org.orangepalantir.genericpca;

import lightgraph.DataSet;
import lightgraph.Graph;
import lightgraph.GraphPoints;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import java.awt.Dimension;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * Created by msmith on 04.10.17.
 */
public class CoefficientKmeansND {
    List<List<IndexedCoefficient>> coefficients;
    int ks = 4;
    int levels = 100;
    List<Path> labels;

    public void setInput(List<List<IndexedCoefficient>> input){
        List<List<IndexedCoefficient>> replacement = new ArrayList<>(input.size());
        for(List<IndexedCoefficient> shape: input){
            List<IndexedCoefficient> sorted = new ArrayList<>(shape);
            sorted.sort((a,b)->Integer.compare(-a.i, -b.i));
            replacement.add(sorted);
        }
        coefficients = replacement;
    }

    static double difference(double[] a, double[] b){
        double sum = 0;
        for(int i = 0; i<a.length; i++){
            sum += (a[i] - b[i])*(a[i] - b[i]);
        }
        return Math.sqrt(sum);
    }

    public void plot(int[] indexes) throws IOException {
        int n = indexes.length;
        double[] data = new double[n*coefficients.size()];
        for(int k = 0; k<coefficients.size(); k++){
            List<IndexedCoefficient> shape = coefficients.get(k);
            for(int s = 0; s<n; s++){
                data[n*k + s] = shape.get(indexes[s]).getCoefficient();
            }
        }

        double[] means = getInitialMean(n, data);
        System.out.println("starting");
        double[] umeans = getMeans(n, means, data);
        double diff = difference(means, umeans);
        means = umeans;
        for(int k = 0; k<levels; k++){
            umeans = getMeans(n, means, data);
            double delta = difference(means, umeans);
            System.out.println(delta);
            means = umeans;
            if(delta/diff < 1e-8) break;
        }

        List<List<double[]>> partitions = new ArrayList<>();
        List<Map<Path, AtomicInteger>> partyLabels = new ArrayList<>();

        for(int i = 0; i<ks; i++){
            partitions.add(new ArrayList<>());
            if(labels!=null){
                partyLabels.add(new HashMap<>());
            }
        }
        int coefficientIndex = 0;

        for(List<IndexedCoefficient> shape: coefficients){
            double[] vector = new double[indexes.length];
            for(int i = 0; i<n; i++){
                vector[i] = shape.get(i).getCoefficient();
            }

            Double min = Double.MAX_VALUE;

            int dex = 0;
            for(int s = 0; s<ks; s++){
                double d = 0;
                for(int i = 0; i<n; i++){
                    double v = (vector[i] - means[s*n + i]);
                    d += v*v;
                }

                if(d<min){
                    dex = s;
                    min = d;
                }

            }

            partitions.get(dex).add(vector);
            if(labels!=null){
                partyLabels.get(dex).computeIfAbsent(labels.get(coefficientIndex).getParent(), a->new AtomicInteger(0)).getAndIncrement();
            }
            coefficientIndex++;
        }

        Graph graph = new Graph();
        int count = 0;
        List<GraphPoints> gp = GraphPoints.getGraphPoints();
        for(List<double[]> part: partitions){
            if(part.size()==0) continue;
            double[] x = new double[part.size()];
            double[] y = new double[part.size()];
            for(int o = 0; o<part.size(); o++){
                x[o] = part.get(o)[0];
                y[o] = part.get(o)[1];
            }
            DataSet set = graph.addData(x, y);
            set.setLine(null);
            set.setPoints(gp.get(count%9 + 3));
            set.setLabel(String.format("x dex: %d, y dex: %d", count/ks, count%ks));
            count++;


        }

        StringBuilder builds = new StringBuilder("Separated on indexes");
        for(int index: indexes){
            builds.append(String.format(" %d", index));
        }

        graph.setTitle(builds.toString());
        graph.show(true);
        if(labels!=null){
            showLabels(indexes, partyLabels);
        }
    }

    /**
     * Gets an auto-generated means for the provided coefficient. Assuming coefficients can be positive
     * or negative
     * @param
     * @return
     */
    double[] getInitialMean(int n, double[] values){
        double[] mins = new double[n];
        double[] maxs = new double[n];

        for(int i = 0; i<n; i++){
            mins[i] = Double.MAX_VALUE;
            maxs[i] = -Double.MAX_VALUE;
        }

        int vectors = values.length/n;

        for(int i = 0; i<vectors; i++){
            for(int j = 0; j<n; j++){
                double v = values[i*n + j];
                mins[j] = v<mins[j]?v:mins[j];
                maxs[j] = v>maxs[j]?v:maxs[j];
            }
        }

        double[] delta = new double[n];
        for(int i = 0; i<n; i++){
            delta[i] = maxs[i] - mins[i];
        }

        double[] initial = new double[ks*n];
        for(int k = 0; k<ks; k++){

            for(int i = 0; i<n; i++){
                initial[k*n + i] = mins[i] + Math.random()*delta[i];
            }

        }


        return initial;
    }

    double[] getMeans(int n, double[] oldmeans, double[] data){

        double[] updated = new double[ks*n];
        double[] counts = new double[ks];
        double[] vector = new double[n];
        int vectors = data.length/n;
        for(int i = 0; i<vectors; i++){
            for(int k = 0; k<n; k++){
                vector[k] = data[n*i + k];
            }

            double min = Double.MAX_VALUE;
            int dex = -1;
            for(int j = 0; j<ks; j++){
                double d = 0;
                for(int k = 0; k<n; k++){
                    double delta = vector[k] - oldmeans[n*j + k];
                    d += delta*delta;
                }

                if(d<min){
                    dex = j;
                    min = d;
                }

            }
            for(int j = 0; j<n; j++){
                updated[dex*n + j] += vector[j];
            }
            counts[dex]++;
        }

        for(int i = 0; i<ks; i++){
            double c = counts[i];
            if(c>0){
                for(int j = 0; j<n; j++){
                    updated[i*n + j] /= c;
                }
            }
        }
        return updated;
    }

    public void showLabels(int[] indexes, List<Map<Path, AtomicInteger>> labelParty) throws IOException {
        String title = String.format("Indexes: %s", Arrays.toString(indexes));
        JFrame frame = new JFrame(title);

        StringBuffer buffer = new StringBuffer();

        int max = 0;
        List<List<Map.Entry<Path, AtomicInteger>>> party = new ArrayList<>();
        for(Map<Path, AtomicInteger> labels: labelParty){
            ArrayList<Map.Entry<Path, AtomicInteger>> entries = new ArrayList<>();
            for(Map.Entry<Path, AtomicInteger> entry: labels.entrySet()){
                entries.add(entry);
            }

            party.add(entries);
            int s = labels.size();
            max = s>max?s:max;
        }

        for(int i = 0; i<max; i++){

            for(List<Map.Entry<Path, AtomicInteger>> labels: party){

                if(i<labels.size()){
                    Map.Entry<Path, AtomicInteger> entry = labels.get(i);
                    buffer.append(entry.getValue() + ":" + entry.getKey());
                }
                buffer.append("\t");

            }
            buffer.append("\n");
        }

        JTextArea area = new JTextArea(buffer.toString());

        frame.add(new JScrollPane(area));

        frame.setSize(new Dimension(800, 800));
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setVisible(true);
        List<String> wrap = new ArrayList<>();
        wrap.add(buffer.toString());
        StringBuilder builds = new StringBuilder("nd_");
        for(int i: indexes){
            builds.append(String.format("-%d", i));
        }
        Files.write(Paths.get(builds.toString()), wrap, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);

    }



    public void setLabels(List<String> labels){
        this.labels = labels.stream().map(Paths::get).collect(Collectors.toList());

    }

    public static void main(String[] args) throws IOException {
        List<List<IndexedCoefficient>> coefficients = IndexedCoefficient.readCoefficients(Paths.get(args[0]));
        CoefficientKmeansND kmeans = new CoefficientKmeansND();
        kmeans.setInput(coefficients);

        if(args.length>=2) {
            List<String> labels = Files.readAllLines(Paths.get(args[1]));
            kmeans.setLabels(labels);
        }

        int x = 1020;
        int y = 1019;
        int z = 1018;
        for(int i = 8; i<9; i++){
            kmeans.ks = i;
            kmeans.plot(new int[]{x, y, z, 1022, 1023, 1021});
            kmeans.plot(new int[]{y, z, x, 1022, 1023, 1021});
            kmeans.plot(new int[]{x, z, y, 1022, 1023, 1021});
        }


    }
}
