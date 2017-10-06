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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * Created by msmith on 04.10.17.
 */
public class CoefficientKmeans {
    List<List<IndexedCoefficient>> coefficients;
    int ks = 4;
    int levels = 10;
    List<Path> labels;

    public void setInput(List<List<IndexedCoefficient>> input){
        List<List<IndexedCoefficient>> replacement = new ArrayList<>(input.size());
        for(List<IndexedCoefficient> shape: input){
            List<IndexedCoefficient> sorted = new ArrayList<>(shape);
            sorted.sort((a,b)->Integer.compare(a.i, b.i));
            replacement.add(sorted);
        }
        coefficients = replacement;
    }

    public void plot(int i, int j) throws IOException {
        double[] ibounds = getBounds(i);
        double[] jbounds = getBounds(j);



        List<List<double[]>> partitions = new ArrayList<>();
        List<Map<Path, AtomicInteger>> partyLabels = new ArrayList<>();

        for(int n = 0; n<ks; n++){
            for(int m = 0; m<ks; m++){
                partitions.add(new ArrayList<>());
                if(labels!=null){
                    partyLabels.add(new HashMap<>());
                }
            }
        }

        int dex = 0;
        for(List<IndexedCoefficient> shape: coefficients){

            int n,m;


            double v = shape.get(i).getCoefficient();
            for(n = 0; n<ks-1; n++){
                if(v<ibounds[n]){
                    break;
                }
            }

            double u = shape.get(j).getCoefficient();
            for(m = 0; m<ks-1; m++){
                if(u<jbounds[m]){
                    break;
                }
            }
            partitions.get(n*ks + m).add(new double[]{v, u});
            if(labels!=null){
                partyLabels.get(n*ks + m).computeIfAbsent(labels.get(dex).getParent(), a->new AtomicInteger(0)).getAndIncrement();
            }
            dex++;
        }

        Graph graph = new Graph();
        int count = 0;
        List<GraphPoints> gp = GraphPoints.getGraphPoints();
        for(List<double[]> part: partitions){
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
        graph.setTitle(String.format("Separated on index %d and index %d", i, j));
        graph.show(true);
        if(labels!=null){
            showLabels(i, j, partyLabels);
        }
    }

    /**
     * Gets an auto-generated boundary for the provided coefficient.
     * @param i
     * @return
     */
    double[] getBounds(int i){
        double[] values = new double[coefficients.size()];
        double min = Double.MAX_VALUE;
        double max = -min;
        int j = 0;
        for(List<IndexedCoefficient> shape: coefficients){
            double v = shape.get(i).getCoefficient();
            values[j] = v;
            min = v<min?v:min;
            max = v>max?v:max;
        }
        double[] boundary = new double[ks-1];
        double delta = (max - min )/ks;

        for(j = 0; j<ks-1; j++){
            boundary[j] = min + (j+1)*delta;
        }

        for(j = 0; j<levels; j++){
            double[] means = getMeans(boundary, values);
            boundary = getBoundary(means);
        }

        return boundary;
    }

    double[] getMeans(double[] boundary, double[] data){

        double[] updated = new double[ks];
        double[] counts = new double[ks];
        pixel:
        for(int i = 0; i<data.length; i++){
            double v = data[i];
            for(int j = 0; j<boundary.length; j++){
                if(v<=boundary[j]){
                    updated[j]+=v;
                    counts[j]++;
                    continue pixel;
                }
            }
            updated[ks-1]+=v;
            counts[ks-1]++;
        }

        for(int i = 0; i<ks; i++){
            double c = counts[i];
            if(c>0){
                updated[i] /= c;
            }
        }
        return updated;
    }

    public void showLabels(int ii, int j, List<Map<Path, AtomicInteger>> labelParty) throws IOException {
        String title = String.format("Separated on index %d and index %d", ii, j);
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
        Files.write(Paths.get(String.format("i%d-j%d.txt", ii, j)), wrap, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);

    }

    double[] getBoundary(double[] means){
        double[] bounds = new double[means.length-1];
        for(int i = 0; i<bounds.length; i++){
            bounds[i] = 0.5*(means[i] + means[i+1]);
        }
        return bounds;
    }

    public void setLabels(List<String> labels){
        this.labels = labels.stream().map(Paths::get).collect(Collectors.toList());

    }

    public static void main(String[] args) throws IOException {
        List<List<IndexedCoefficient>> coefficients = IndexedCoefficient.readCoefficients(Paths.get(args[0]));
        CoefficientKmeans kmeans = new CoefficientKmeans();
        kmeans.setInput(coefficients);
        if(args.length>=2) {
            List<String> labels = Files.readAllLines(Paths.get(args[1]));
            kmeans.setLabels(labels);
        }
        kmeans.plot(1020, 1020);
        kmeans.plot(1019, 1019);
        kmeans.plot(1015, 1015);
        kmeans.plot(1015, 1019);
    }
}
