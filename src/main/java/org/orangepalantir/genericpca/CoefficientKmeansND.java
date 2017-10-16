package org.orangepalantir.genericpca;

import lightgraph.DataSet;
import lightgraph.Graph;
import lightgraph.GraphPoints;
import org.orangepalantir.genericpca.display.TwoDHeatMap;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Created by msmith on 04.10.17.
 */
public class CoefficientKmeansND {
    List<List<IndexedCoefficient>> coefficients;
    int ks = 4;
    int levels = 1000;
    List<Path> labels;
    List<double[]> eigenVectors;
    List<Highlight> highlights;

    public void setInput(List<List<IndexedCoefficient>> input){
        List<List<IndexedCoefficient>> replacement = new ArrayList<>(input.size());
        for(List<IndexedCoefficient> shape: input){
            List<IndexedCoefficient> sorted = new ArrayList<>(shape);
            sorted.sort((a,b)->Integer.compare(a.i, b.i));
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

    public double plot(int[][] groups) throws IOException {
        int[] indexes = new int[groups.length];
        for(int i = 0; i<groups.length; i++){
            indexes[i] = groups[i][0];
        }
        int n = indexes.length;
        double[] data = new double[n*coefficients.size()];
        for(int k = 0; k<coefficients.size(); k++){
            List<IndexedCoefficient> shape = coefficients.get(k);
            for(int s = 0; s<n; s++){
                int[] group = groups[s];
                for(int l = 0; l<group.length; l++){
                    IndexedCoefficient ic = shape.get(group[l]);
                    if(ic.i!=group[l]) throw new RuntimeException("coefficient does not correspond to index!");
                    data[n*k + s] += shape.get(group[l]).getCoefficient();
                }
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
        List<Map<Path, List<double[]>>> partyLabels = new ArrayList<>();

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
                vector[i] = shape.get(indexes[i]).getCoefficient();
                int[] group = groups[i];
                for(int l = 0; l<group.length; l++){
                    vector[i] += shape.get(group[l]).getCoefficient();
                }
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
                partyLabels.get(dex).computeIfAbsent(labels.get(coefficientIndex).getParent(), a->new ArrayList<>()).add(vector);
            }
            coefficientIndex++;
        }

        Graph graph = new Graph();
        int count = 0;

        List<GraphPoints> gp = GraphPoints.getGraphPoints();
        for(List<double[]> part: partitions){
            if(part.size()==0){
                count++;
                continue;
            }
            double[] x = new double[part.size()];
            double[] y = new double[part.size()];
            for(int o = 0; o<part.size(); o++){
                x[o] = part.get(o)[0];
                y[o] = part.get(o)[1];
            }
            DataSet set = graph.addData(x, y);
            set.setLine(null);
            set.setPoints(gp.get(count%9 + 3));
            set.setLabel(String.format("k: %d", count));
            count++;


        }
        if(highlights!=null){
            Color[] colors = {
                    Color.BLUE,
                    Color.RED,
                    Color.ORANGE,
                    Color.DARK_GRAY
            };
            int counter = 0;
            for(Highlight high: highlights){
                List<double[]> values = new ArrayList<>();
                for(Path p: high.conditions){
                    for(Map<Path, List<double[]>> map: partyLabels){
                        values.addAll(map.getOrDefault(p, Collections.EMPTY_LIST));
                    }
                }
                double[] x = new double[values.size()];
                double[] y = new double[values.size()];
                for(int i = 0; i<x.length; i++){
                    x[i] = values.get(i)[0];
                    y[i] = values.get(i)[1];
                }
                DataSet set = graph.addData(x, y);
                set.setLine(null);
                set.setPoints(GraphPoints.crossX());
                set.setColor(colors[counter]);
                set.setPointWeight(2.0);
                counter = (counter+1)%colors.length;
                set.setLabel(high.label);
            }
        }

        StringBuilder builds = new StringBuilder("Separated on indexes");
        for(int index: indexes){
            builds.append(String.format(" %d", index));
        }

        graph.setTitle(builds.toString());
        graph.show(false);
        if(labels!=null){
            showLabels(indexes, partyLabels, means);
        }

        if(eigenVectors!=null){
            writeMeanShapes(indexes, means);
        }

        return calculateVariation(means, partitions, n);
    }

    private double calculateVariation(double[] means, List<List<double[]>> partitions, int n) {
        double sigma = 0;
        for(int i = 0; i<ks; i++){
            List<double[]> vectors = partitions.get(i);
            for(double[] vector: vectors){

                for(int j = 0; j<n; j++){
                    double delta = vector[j] - means[i*n + j];
                    sigma += delta*delta;
                }

            }
        }
        return sigma;
    }

    /**
     * Creates a mean shape based on the clustered coefficients.
     *
     * @param indexes
     * @param means
     * @param coefficients
     */
    void writeMeanShapes(int[] indexes, double[] means){
        int space= eigenVectors.get(0).length;
        int width = (int)Math.sqrt(space);


        for(int j = 0; j<ks; j++){
            StringBuilder name = new StringBuilder("km");
            name.append(j);
            name.append("_");
            for(int dex: indexes){
                name.append(dex);
                name.append("-");
            }
            name.append("sum.png");
            int n = indexes.length;
            double[] output = new double[space];
            for(int i = 0; i<n; i++){
                double[] ev = eigenVectors.get(indexes[i]);
                double ai = means[j*n + i];
                for(int k = 0; k<space; k++){
                    output[k] += ev[k]*ai;
                }
            }
            BufferedImage img = TwoDHeatMap.createMap(width, width, 10, 10, output);

            try {
                ImageIO.write(img, "PNG", new File(name.toString()));
            } catch (IOException e) {
                e.printStackTrace();
            }


        }
    }

    /**
     * Partition by means.
     *
     */
    public List<List<List<IndexedCoefficient>>> partitionByMeans(int[] indexes, double[] means, List<List<IndexedCoefficient>> coefficients){
        return null;
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
    static class Highlight{
        String label;
        List<Path> conditions;
        public Highlight(String all){
            String[] tokens = all.split(Pattern.quote("\t"));
            label = tokens[0];
            conditions = new ArrayList<>();
            for(int i = 1; i<tokens.length; i++){
                conditions.add(Paths.get(tokens[i]));
            }
        }
    }
    public void setHighlights(List<String> highlights){
        System.out.println("working");
        this.highlights = highlights.stream().map(Highlight::new).collect(Collectors.toList());
    }

    public void showLabels(int[] indexes, List<Map<Path, List<double[]>>> labelParty, double[] means) throws IOException {
        String title = String.format("Indexes: %s", Arrays.toString(indexes));
        JFrame frame = new JFrame(title);

        StringBuffer buffer = new StringBuffer();

        int max = 0;
        int n = indexes.length;
        for(int i = 0; i<ks; i++){
            buffer.append("#k" + i);
            for(int j = 0; j<n; j++){
                buffer.append(String.format("\t%f", means[i*n + j]));
            }
            buffer.append("\n");
        }
        List<List<Map.Entry<Path, List<double[]>>>> party = new ArrayList<>();
        for(Map<Path, List<double[]>> labels: labelParty){
            ArrayList<Map.Entry<Path, List<double[]>>> entries = new ArrayList<>();
            for(Map.Entry<Path, List<double[]>> entry: labels.entrySet()){
                entries.add(entry);
            }

            party.add(entries);
            int s = labels.size();
            max = s>max?s:max;
        }

        for(int i = 0; i<max; i++){

            for(List<Map.Entry<Path, List<double[]>>> labels: party){

                if(i<labels.size()){
                    Map.Entry<Path, List<double[]>> entry = labels.get(i);
                    buffer.append(entry.getValue().size() + ":" + entry.getKey());
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
        builds.append(".txt");
        Files.write(Paths.get(builds.toString()), wrap, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);

    }
    public void setEigens(List<String> eigens){
        eigenVectors = new ArrayList<double[]>();
        for(String line: eigens){
            double[] values = Arrays.stream(line.split(Pattern.quote("\t"))).mapToDouble(Double::parseDouble).toArray();
            eigenVectors.add(values);
        }
    }


    public void setLabels(List<String> labels){
        this.labels = labels.stream().map(Paths::get).collect(Collectors.toList());

    }

    public static void main(String[] args) throws IOException {
        //testVersion();
        loadDataAndRun(args);
    }

    public static void loadDataAndRun(String[] args) throws IOException {
        List<List<IndexedCoefficient>> coefficients = IndexedCoefficient.readCoefficients(Paths.get(args[0]));
        CoefficientKmeansND kmeans = new CoefficientKmeansND();
        kmeans.setInput(coefficients);

        if(args.length>=2) {
            List<String> labels = Files.readAllLines(Paths.get(args[1]));
            kmeans.setLabels(labels);
        }
        if(args.length>=3){
            List<String> eigens = Files.readAllLines(Paths.get(args[2]));
            kmeans.setEigens(eigens);
        }

        if(args.length>=4){
            List<String> highlights = Files.readAllLines(Paths.get(args[3]));
            kmeans.setHighlights(highlights);
        }

        double[] x = new double[10];
        double[] y = new double[10];

        for(int i = 0; i<10; i++){
            int top = 255-2*i;
            int[][] indexes = {
                    { top},
                    {top - 1}
            };
            kmeans.ks = 6;
            //241 & 240 ntc separates.
            double s = kmeans.plot(indexes);
            x[i] = kmeans.ks;
            y[i] = s;
        }

        new Graph(x, y).show(true, "variance");
    }

    public static void testVersion(){

        double[] v1 = {5*Math.sqrt(2)/2, 5*Math.sqrt(2)/2, 0};
        double[] v2 = {-3*Math.sqrt(2)/2, 3*Math.sqrt(2)/2, 0};
        double[] v3 = {1, 0, 3};
        List<List<IndexedCoefficient>> data = new ArrayList<>(300);
        List<String> labels = new ArrayList<>(300);
        Random ng = new Random();

        double noise = 1.0;

        for(int i = 0; i<100; i++){
            List<IndexedCoefficient> shape = new ArrayList<>(2);
            for(int j = 0; j<3; j++){
                IndexedCoefficient ic =new IndexedCoefficient(j, v1[j] + ng.nextGaussian()*noise);
                shape.add(ic);
            }
            labels.add("one/d.txt");
            data.add(shape);
        }

        for(int i = 0; i<100; i++){
            List<IndexedCoefficient> shape = new ArrayList<>(2);
            for(int j = 0; j<3; j++){
                IndexedCoefficient ic =new IndexedCoefficient(j, v2[j] + ng.nextGaussian()*noise);
                shape.add(ic);
            }
            labels.add("two/d.txt");
            data.add(shape);
        }

        for(int i = 0; i<100; i++){
            List<IndexedCoefficient> shape = new ArrayList<>(2);
            for(int j = 0; j<3; j++){
                IndexedCoefficient ic =new IndexedCoefficient(j, v3[j] + ng.nextGaussian()*noise);
                shape.add(ic);
            }
            labels.add("three/d.txt");
            data.add(shape);
        }

        CoefficientKmeansND kmeans = new CoefficientKmeansND();
        kmeans.setInput(data);
        kmeans.setLabels(labels);
        double[] x = new double[5];
        double[] y = new double[5];
        try {
            for(int i = 1; i<=5; i++){
                kmeans.ks = i;
                double sigma = kmeans.plot(new int[][]{{0}, {1}, {2}});
                x[i-1] = i;
                y[i-1] = sigma;
            }
            new Graph(x, y).show(true);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
