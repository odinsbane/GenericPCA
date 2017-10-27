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
    List<List<IndexedCoefficient>> original;
    List<List<IndexedCoefficient>> normalized;
    int ks = 4;
    int levels = 1000;
    List<Path> labels;
    List<double[]> eigenVectors;
    List<Highlight> highlights;
    static boolean unscaled = true;
    public void setInput(List<List<IndexedCoefficient>> input){
        List<List<IndexedCoefficient>> replacement = new ArrayList<>(input.size());
        for(List<IndexedCoefficient> shape: input){
            List<IndexedCoefficient> sorted = new ArrayList<>(shape);
            sorted.sort((a,b)->Integer.compare(a.i, b.i));
            replacement.add(sorted);
        }
        original = replacement;
        normalized = normalizeCoefficients(replacement);
    }

    static double difference(double[] a, double[] b){
        double sum = 0;
        for(int i = 0; i<a.length; i++){
            sum += (a[i] - b[i])*(a[i] - b[i]);
        }
        return Math.sqrt(sum);
    }

    public double calculate(int[] indexes){
        int n = indexes.length;
        double[] data = new double[n*normalized.size()];
        for(int k = 0; k<normalized.size(); k++){
            List<IndexedCoefficient> shape = normalized.get(k);
            for(int s = 0; s<n; s++){
                IndexedCoefficient ic = shape.get(indexes[s]);
                if(ic.i!=indexes[s]) throw new RuntimeException("coefficient does not correspond to index!");
                data[n*k + s] += ic.getCoefficient();
            }
        }

        double[] means = getInitialMean(n, data);
        double[] umeans = getMeans(n, means, data);
        double diff = difference(means, umeans);
        means = umeans;
        for(int k = 0; k<levels; k++){
            umeans = getMeans(n, means, data);
            double delta = difference(means, umeans);
            means = umeans;
            if(delta==0) break;
        }
        umeans = getMeans(n, means, data);

        List<List<double[]>> partitions = new ArrayList<>();

        for(int i = 0; i<ks; i++){
            partitions.add(new ArrayList<>());
        }

        for(List<IndexedCoefficient> shape: normalized){
            double[] vector = new double[indexes.length];
            for(int i = 0; i<n; i++){
                vector[i] += shape.get(indexes[i]).getCoefficient();
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
        }
        return calculateVariation(means, partitions, n);
    }
    private List<List<IndexedCoefficient>> normalizeCoefficients(List<List<IndexedCoefficient>> co){
        int cnets = co.get(0).size();
        List<List<IndexedCoefficient>> normed = new ArrayList<>(co.size());
        for(List<IndexedCoefficient> shape: co){
            normed.add(new ArrayList<>());
        }
        for(int index = 0; index< cnets; index++){
            double min = Double.MAX_VALUE;
            double max = -Double.MAX_VALUE;
            double sum = 0;
            double sum2 = 0;
            for(List<IndexedCoefficient> shape: co){

                double v = shape.get(index).getCoefficient();
                sum += v;
                sum2 += v*v;
                min = v<min?v:min;
                max = v>max?v:max;

            }
            sum = sum/co.size();
            sum2 = Math.sqrt(sum2/co.size() - sum*sum);
            min = sum - sum2;
            max = sum + sum2;
            for(int i = 0; i<co.size(); i++){
                List<IndexedCoefficient> shape = co.get(i);
                List<IndexedCoefficient> nShape = normed.get(i);
                IndexedCoefficient ic = shape.get(index);
                double v = ic.getCoefficient();
                double scaled = 2*(v - min)/(max - min) - 1;
                nShape.add(new IndexedCoefficient(ic.i, scaled));
            }

        }
        return normed;
    }

    public double plot(int[] indexes) throws IOException {
        int n = indexes.length;
        double[] data = new double[n*normalized.size()];
        for(int k = 0; k<normalized.size(); k++){
            List<IndexedCoefficient> shape = normalized.get(k);
            for(int s = 0; s<n; s++){
                IndexedCoefficient ic = shape.get(indexes[s]);
                if(ic.i!=indexes[s]) throw new RuntimeException("coefficient does not correspond to index!");
                data[n*k + s] += ic.getCoefficient();
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
            if(delta==0) break;
        }
        umeans = getMeans(n, means, data);

        List<List<double[]>> partitions = new ArrayList<>();
        List<Map<Path, List<double[]>>> partyLabels = new ArrayList<>();

        for(int i = 0; i<ks; i++){
            partitions.add(new ArrayList<>());
            if(labels!=null){
                partyLabels.add(new HashMap<>());
            }
        }
        int coefficientIndex = 0;

        for(List<IndexedCoefficient> shape: normalized){
            double[] vector = new double[indexes.length];
            for(int i = 0; i<n; i++){
                vector[i] += shape.get(indexes[i]).getCoefficient();
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
        double xmin = Double.MAX_VALUE;
        double ymin = Double.MAX_VALUE;
        double xmax = -xmin;
        double ymax = -ymin;
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

            double[] xrange = getAxisExtremes(x);
            double[] yrange = getAxisExtremes(y);
            xmin = xmin<xrange[0]?xmin:xrange[0];
            xmax = xmax>xrange[1]?xmax:xrange[1];
            ymin = ymin<yrange[0]?ymin:yrange[0];
            ymax = ymax>yrange[1]?ymax:yrange[1];

            set.setLine(null);
            //set.setPoints(gp.get(count%9 + 3));
            set.setPoints(GraphPoints.filledCircles());
            Color c = set.COLOR;
            c = new Color(c.getRed(), c.getGreen(), c.getBlue(), 050);
            set.setColor(c);
            set.setLabel(String.format("k: %d", count));
            count++;


        }
        graph.setXRange(xmin, xmax);
        graph.setYRange(ymin, ymax);
        if(highlights!=null){
            Color[] colors = {
                    Color.BLUE,
                    Color.RED,
                    Color.DARK_GRAY
            };
            GraphPoints[] points = {
                    GraphPoints.hollowSquares(),
                    GraphPoints.hollowCircles(),
                    GraphPoints.hollowTriangles()
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
                set.setPoints(points[counter%points.length]);
                set.setColor(colors[(counter/colors.length)%colors.length]);
                set.setPointWeight(2.0);
                set.setPointSize(8.0);
                counter++;
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
        double[] single = new double[n];
        double[] singleSqd = new double[n];
        double counter = 0;
        double massSigma = 0;
        double sigma = 0;
        for(int i = 0; i<ks; i++){
            List<double[]> vectors = partitions.get(i);

            for(double[] vector: vectors){

                for(int j = 0; j<n; j++){
                    double delta = vector[j] - means[i*n + j];
                    sigma += delta*delta;
                    single[j] += vector[j];
                    singleSqd[j] += vector[j]*vector[j];
                }
                counter++;
            }
        }
        sigma = sigma/counter;

        for(int i = 0; i<n; i++){
            double xbar = single[i]/counter;
            massSigma += singleSqd[i]/counter - xbar*xbar;
        }

        return sigma/massSigma;
    }
    /**
     * Creates a mean shape based on the clustered coefficients.
     *
     * @param indexes
     * @param means
     */
    void writeMeanShapes(int[] indexes, double[] means){
        int space= eigenVectors.get(0).length;
        int width = (int)Math.sqrt(space);
        List<List<List<IndexedCoefficient>>> partitionedShapes = new ArrayList<>();
        int n = indexes.length;

        for(int i = 0; i<ks; i++){
            partitionedShapes.add(new ArrayList<>());
        }

        for(int shapeIndex = 0; shapeIndex<normalized.size(); shapeIndex++){
            List<IndexedCoefficient> shape = normalized.get(shapeIndex);
            Double min = Double.MAX_VALUE;

            int dex = 0;
            for(int s = 0; s<ks; s++){
                double d = 0;
                for(int i = 0; i<n; i++){
                    double v = (shape.get(indexes[i]).getCoefficient() - means[s*n + i]);
                    d += v*v;
                }

                if(d<min){
                    dex = s;
                    min = d;
                }

            }

            partitionedShapes.get(dex).add(original.get(shapeIndex));

        }
        for(int j = 0; j<ks; j++){
            StringBuilder name = new StringBuilder("km");
            name.append("_");
            for(int dex: indexes){
                name.append(dex);
                name.append("-");
            }
            name.append(j);
            name.append("sum.png");
            double[] output = new double[space];
            double[] a = new double[eigenVectors.size()];
            List<List<IndexedCoefficient>> partition = partitionedShapes.get(j);
            for(List<IndexedCoefficient> shape: partition){
                for(int i = 0; i<a.length; i++){
                    a[i] += shape.get(i).getCoefficient();
                }
            }
            int count = partition.size();
            for(int i = 0; i<a.length; i++){
                a[i] = a[i]/count;
            }
            for(int i = 0; i<a.length; i++){
                double[] ev = eigenVectors.get(i);
                double ai = a[i];
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
        Random rand = new Random();
        for(int k = 0; k<ks; k++){

            for(int i = 0; i<n; i++){
                initial[k*n + i] = mins[i] + rand.nextDouble()*delta[i];
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
    static public double[] getAxisExtremes(double[] data){
        if(unscaled) return new double[] {-3, 3};
        double sum = 0;
        double sum2 = 0;
        double min = Double.MAX_VALUE;
        double max = -min;
        for(double d: data){
            sum += d;
            sum2 += d*d;
            min = d<min?d:min;
            max = d>max?d:max;


        }

        double center = sum/data.length;
        double stdev = Math.sqrt(sum2/data.length - center*center);
        double top = center + 2*stdev;
        top = top>max?max:top;
        double bottom = center - 2*stdev;
        bottom = bottom<min?min:bottom;

        return new double[]{bottom, top};


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
        kmeans.normalizeCoefficients(coefficients);
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


        int ks = 1;
        Graph variancePlot = new Graph();
        int[][] indexGroups = new int[][]{
                {1023, 1022},
                {1021, 1020},
                {1019, 1018},
                {1017, 1016},
                {1015, 1014},
                {1013, 1012}
        };
        int vectors = indexGroups.length;

        for(int i = 0; i<vectors; i++){
            int prime = 66;
            int aux =  171;
            int[] indexes = indexGroups[i];

            double[] x = new double[ks];
            double[] y = new double[ks];

            for(int j = 0; j<ks; j++) {
                kmeans.ks = 6 + j*3;
                //241 & 240 ntc separates.
                double s = kmeans.plot(indexes);
                x[j] = kmeans.ks;
                y[j] = s;
            }
            DataSet set = variancePlot.addData(x, y);
            set.setLabel(Arrays.toString(indexes));
        }
        variancePlot.show(true, "variance");

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
                double sigma = kmeans.plot(new int[]{0, 1, 2});
                x[i-1] = i;
                y[i-1] = sigma;
            }
            new Graph(x, y).show(true);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
