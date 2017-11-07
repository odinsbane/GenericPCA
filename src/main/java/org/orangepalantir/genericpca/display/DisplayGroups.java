package org.orangepalantir.genericpca.display;

import lightgraph.DataSet;
import lightgraph.Graph;
import lightgraph.GraphDefaults;
import lightgraph.GraphPoints;
import lightgraph.painters.GraphPainter;
import org.orangepalantir.genericpca.IndexedCoefficient;
import org.orangepalantir.genericpca.SimpleReader;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import java.awt.Color;
import java.awt.EventQueue;
import java.awt.FileDialog;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Created by msmith on 02.11.17.
 */
public class DisplayGroups {
    Graph display;
    Map<String, DisplayGroup> groups = new HashMap<>();
    List<double[]> data;
    List<String> labels;
    Map<String, Set<String>> highlights;
    int xDex = -1;
    int yDex = -1;
    public static String DEFAULT_KEY = "default";
    public void buildGui(){

        Path coefficients, rubric, highlight;
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.PAGE_AXIS));
        JLabel coefficientLabel = new JLabel("-none-");
        JLabel labelLabel = new JLabel("-none-");
        JLabel groupLabel = new JLabel("-none-");
        JTextField xPlotDex = new JTextField("-1");
        JTextField yPlotDex = new JTextField("-1");


        panel.add(coefficientLabel);
        panel.add(labelLabel);
        panel.add(groupLabel);
        panel.add(new JLabel("first index"));
        panel.add(xPlotDex);
        panel.add(new JLabel("second index"));
        panel.add(yPlotDex);

        JButton process = new JButton("process");

        JFrame frame = new JFrame("selector");

        coefficientLabel.addMouseListener(new MouseAdapter(){
            @Override
            public void mouseClicked(MouseEvent evt){
                FileDialog fd = new FileDialog(frame, "Select file with coefficient data");
                fd.setMode(FileDialog.LOAD);

                fd.setVisible(true);
                String file = fd.getFile();
                String dir = fd.getDirectory();
                if(file!=null){
                    try {
                        loadData(Paths.get(dir, file));
                        int max = data.get(0).length - 1;
                        int next = max - 1;
                        xDex = max;
                        yDex = next;
                        xPlotDex.setText("" + xDex);
                        yPlotDex.setText("" + yDex);

                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    coefficientLabel.setText(file);
                }
            }
        });

        groupLabel.addMouseListener(new MouseAdapter(){
            @Override
            public void mouseClicked(MouseEvent evt){
                FileDialog fd = new FileDialog(frame, "Select file with highlight information data");
                fd.setMode(FileDialog.LOAD);

                fd.setVisible(true);
                String file = fd.getFile();
                String dir = fd.getDirectory();
                if(file!=null){
                    try {
                        loadHighlights(Paths.get(dir, file));
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    groupLabel.setText(file);
                }
            }
        });
        labelLabel.addMouseListener(new MouseAdapter(){
            @Override
            public void mouseClicked(MouseEvent evt){
                FileDialog fd = new FileDialog(frame, "Select file with highlight information data");
                fd.setMode(FileDialog.LOAD);

                fd.setVisible(true);
                String file = fd.getFile();
                String dir = fd.getDirectory();
                if(file!=null){
                    try {
                        loadLabels(Paths.get(dir, file));
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    labelLabel.setText(file);
                }
            }
        });

        panel.add(process);

        process.addActionListener(evt->{
            buildData();
            String[] values = groups.keySet().toArray(new String[groups.keySet().size()]);
            JList<String> list = new JList<>(values);
            list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
            list.addListSelectionListener(lll->{
                List<String> selected = list.getSelectedValuesList();
                for(DisplayGroup group: groups.values()){
                    group.setHighlighted(false);
                }
                for(String sel: selected){
                    groups.get(sel).setHighlighted(true);
                }
                display.refresh(false);
            });
            panel.add(list);

        });
        xPlotDex.addActionListener(evt->{
            int i = Integer.parseInt(xPlotDex.getText());
            xDex = i;
        });
        yPlotDex.addActionListener(evt->{
            int i = Integer.parseInt(xPlotDex.getText());
            yDex = i;
        });

        frame.setContentPane(panel);
        frame.pack();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    }

    public void loadData(Path p) throws IOException {
        List<List<IndexedCoefficient>> shapes = IndexedCoefficient.readCoefficients(p);
        data = shapes.stream().map(
                l->l.stream().mapToDouble(IndexedCoefficient::getCoefficient).toArray()
        ).collect(Collectors.toList());
    }

    public void loadLabels(Path p) throws IOException {
        labels = Files.readAllLines(p);
    }

    public void loadHighlights(Path p) throws IOException {
        highlights = new HashMap<>();
        try(BufferedReader reader = Files.newBufferedReader(p, StandardCharsets.UTF_8)){
            String line;
            while((line=reader.readLine())!=null) {
                String[] tokens = line.split("\\t");
                highlights.put(tokens[0], new HashSet<>());
                for (int i = 1; i < tokens.length; i++) {
                    highlights.get(tokens[0]).add(tokens[i]);
                }
            }
        }
    }

    class DisplayGroup{
        String label;
        List<double[]> dataPoints = new ArrayList<>();
        List<String> labels = new ArrayList<>();
        double[] x, y;
        boolean highlighted = false;
        public DisplayGroup(String l){
            label = l;
        }

        public void prepareData(){
            x = new double[dataPoints.size()];
            y = new double[dataPoints.size()];
            int i = 0;
            for(double[] row: dataPoints){
                x[i] = row[xDex];
                y[i] = row[yDex];
                i++;
            }
        }
        public double[] getX(){
            return x;
        }

        public double[] getY(){
            return y;
        }

        public void addPoint(String label, double[] values) {
            labels.add(label);
            dataPoints.add(values);
        }

        public boolean highlighted() {
            return highlighted;
        }

        public void setHighlighted(boolean v){
            highlighted = v;
        }
    }

    public List<DisplayGroup> getGroups(){
        groups = new HashMap<>();
        for(String key: highlights.keySet()){
            groups.put(key, new DisplayGroup(key));
        }
        DisplayGroup def = new DisplayGroup(DEFAULT_KEY);

        for(int i = 0; i<data.size(); i++){
            double[] values = data.get(i);
            String label = labels.get(i);
            String condition = getCondition(label);
            boolean grouped = false;
            for(Map.Entry<String, Set<String>> entry: highlights.entrySet()){
                if(entry.getValue().contains(condition)){
                    groups.get(entry.getKey()).addPoint(label, values);
                    grouped = true;
                    break;
                }

            }
            if(!grouped) {
                def.addPoint(label, values);
            }
        }
        groups.put(DEFAULT_KEY, def);
        return new ArrayList<>(groups.values());
    }

    public String getCondition(String s) {
        return s.substring(0, s.lastIndexOf('/'));
    }

    public void buildData(){
        //everything ready?
        if(data==null || labels==null || highlights == null){
            return;
        }

        if(xDex<0 || yDex<0){
            return;
        }
        if(display!=null){
            //data is built and changes were made.
            updateData();
        }else{

            display = new Graph();
            List<DisplayGroup> groups = getGroups();
            for(final DisplayGroup group: groups){
                group.prepareData();
                DataSet set = display.addData(group.getX(), group.getY());
                set.setLabel(set.label);
                set.setLine(null);
                Color c = GraphDefaults.getDefaultColor(display.dataSetCount());
                set.setColor(new Color(c.getRed(), c.getGreen(), c.getBlue(), 50));
                set.setPoints(new GraphPoints() {
                    GraphPoints normal = GraphPoints.hollowCircles();
                    GraphPoints highlight = GraphPoints.filledCircles();
                    {
                        setSize(5);
                    }

                    @Override
                    public void drawPoint(Point2D pt, GraphPainter painter) {

                        if(group.highlighted()){
                            Color c = painter.getColor();
                            Color n = new Color(c.getRed(), c.getGreen(), c.getBlue()).darker();
                            painter.setColor(n);
                            highlight.drawPoint(pt, painter);
                            painter.setColor(c);
                        }

                        normal.drawPoint(pt, painter);
                    }

                    @Override
                    public void setSize(double s){
                        System.out.println("sized" + s);
                        normal.setSize(s);
                        highlight.setSize(2*s);
                    }
                });
            }

            display.show(true, "Scatter Plot " + xDex + " vs " + yDex);
        }

    }

    public void updateData(){

    }
    public static void main(String[] args){
        EventQueue.invokeLater(()->new DisplayGroups().buildGui());
    }

}
