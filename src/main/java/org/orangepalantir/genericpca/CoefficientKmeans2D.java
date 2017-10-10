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
public class CoefficientKmeans2D {
    List<List<IndexedCoefficient>> coefficients;
    int ks = 4;
    int levels = 100;
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

    static double difference(double[] a, double[] b){
        double sum = 0;
        for(int i = 0; i<a.length; i++){
            sum += (a[i] - b[i])*(a[i] - b[i]);
        }
        return Math.sqrt(sum);
    }

    public void plot(int i, int j) throws IOException {

        double[] data = new double[2*coefficients.size()];
        for(int k = 0; k<coefficients.size(); k++){
            data[2*k] = coefficients.get(k).get(i).getCoefficient();
            data[2*k+1] = coefficients.get(k).get(j).getCoefficient();
        }

        double[] means = getInitialMean(i, j);
        double[] umeans = getMeans(means, data);
        double diff = difference(means, umeans);
        means = umeans;
        for(int k = 0; k<levels; k++){
            umeans = getMeans(means, data);
            double delta = difference(means, umeans);
            means = umeans;
            if(delta/diff < 1e-8) break;
        }

        List<List<double[]>> partitions = new ArrayList<>();
        List<Map<Path, AtomicInteger>> partyLabels = new ArrayList<>();

        for(int n = 0; n<ks; n++){
            partitions.add(new ArrayList<>());
            if(labels!=null){
                partyLabels.add(new HashMap<>());
            }
        }
        int coefficientIndex = 0;
        for(List<IndexedCoefficient> shape: coefficients){

            double vx = shape.get(i).getCoefficient();
            double vy = shape.get(j).getCoefficient();
            Double min = Double.MAX_VALUE;
            int dex = 0;
            for(int s = 0; s<ks; s++){
                double mx = means[2*s];
                double my = means[2*s + 1];

                double d = (vx - mx)*(vx -mx) + (vy - my)*(vy - my);

                if(d<min){
                    dex = s;
                    min = d;
                }

            }

            partitions.get(dex).add(new double[]{vx, vy});
            if(labels!=null){
                partyLabels.get(dex).computeIfAbsent(labels.get(coefficientIndex).getParent(), a->new AtomicInteger(0)).getAndIncrement();
            }
            coefficientIndex++;
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
        graph.setTitle(String.format("Separated on index %d and index %d; levels: %d", i, j, levels));
        graph.show(true);
        if(labels!=null){
            showLabels(i, j, partyLabels);
        }
    }

    /**
     * Gets an auto-generated means for the provided coefficient. Assuming coefficients can be positive
     * or negative
     * @param i, j
     * @return
     */
    double[] getInitialMean(int i, int j){
        double minx = Double.MAX_VALUE;
        double miny = Double.MAX_VALUE;

        double maxx = -minx;
        double maxy = -miny;
        int k  = 0;
        for(List<IndexedCoefficient> shape: coefficients){
            double vx = shape.get(i).getCoefficient();
            double vy = shape.get(j).getCoefficient();
            minx = vx<minx?vx:minx;
            maxx = vx>maxx?vx:maxx;
            miny = vy<miny?vy:miny;
            maxy = vy>maxy?vy:maxy;
            k++;
        }

        double deltax = (maxx - minx );
        double deltay = (maxy - miny);

        double dtheta = Math.PI*2/ks;
        double theta0 = deltax>deltay?Math.PI*0.5:0;
        double[] initial = new double[ks*2];
        for(k = 0; k<ks; k++){
            double sin = Math.sin(dtheta*k + theta0);
            double cos = Math.sin(dtheta*k + theta0);
            initial[k*2] = minx + deltax*sin;
            initial[k*2+1] = miny + deltay*cos;

        }


        return initial;
    }

    double[] getMeans(double[] oldmeans, double[] data){

        double[] updated = new double[ks*2];
        double[] counts = new double[ks];
        for(int i = 0; i<data.length/2; i++){
            double vx = data[2*i];
            double vy = data[2*i + 1];

            double min = Double.MAX_VALUE;
            int dex = -1;
            for(int j = 0; j<ks; j++){
                double mx = oldmeans[2*j];
                double my = oldmeans[2*j + 1];

                double d = (vx - mx)*(vx -mx) + (vy - my)*(vy - my);

                if(d<min){
                    dex = j;
                    min = d;
                }

            }
            updated[dex*2] += vx;
            updated[dex*2 + 1] += vy;
            counts[dex]++;
        }

        for(int i = 0; i<ks; i++){
            double c = counts[i];
            if(c>0){
                updated[2*i] /= c;
                updated[2*i+1] /= c;
            }
        }
        return updated;
    }

    public void showLabels(int ii, int j, List<Map<Path, AtomicInteger>> labelParty) throws IOException {
        String title = String.format("2D index %d and index %d levels=%d", ii, j, levels);
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
        Files.write(Paths.get(String.format("d2_i%d-j%d.txt", ii, j)), wrap, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);

    }



    public void setLabels(List<String> labels){
        this.labels = labels.stream().map(Paths::get).collect(Collectors.toList());

    }

    public static void main(String[] args) throws IOException {
        List<List<IndexedCoefficient>> coefficients = IndexedCoefficient.readCoefficients(Paths.get(args[0]));
        CoefficientKmeans2D kmeans = new CoefficientKmeans2D();
        kmeans.setInput(coefficients);

        if(args.length>=2) {
            List<String> labels = Files.readAllLines(Paths.get(args[1]));
            kmeans.setLabels(labels);
        }

        for(int i = 2; i<8; i++){
            kmeans.ks = i;
            kmeans.plot(1015, 1015);
        }


    }
}
