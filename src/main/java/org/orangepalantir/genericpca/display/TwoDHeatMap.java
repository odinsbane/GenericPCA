package org.orangepalantir.genericpca.display;

import org.orangepalantir.genericpca.SimpleReader;
import org.orangepalantir.genericpca.Trainer;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JSlider;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.EventQueue;
import java.awt.FileDialog;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.SinglePixelPackedSampleModel;
import java.awt.image.WritableRaster;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by msmith on 03.08.17.
 */
public class TwoDHeatMap {
    int width;
    int height;
    int xsize = 2;
    int ysize = 2;
    public TwoDHeatMap(int width, int height){
        this.width = width;
        this.height = height;
    }
    public void displayData(List<double[]> values){

        JFrame frame = new JFrame("components");

        JLabel label = new JLabel("waiting");
        JSlider slider = new JSlider();
        slider.setMaximum(values.size()-1);
        slider.setMinimum(0);

        JLabel eigenLabel = new JLabel("waiting");
        JSlider eigenSlider = new JSlider();
        eigenSlider.setMaximum(values.get(0).length-1);
        eigenSlider.setMinimum(0);

        JPanel c = new JPanel();
        c.setLayout(new BorderLayout());
        c.add(label, BorderLayout.CENTER);
        c.add(slider, BorderLayout.SOUTH);

        JPanel c2 = new JPanel();
        c2.setLayout(new BorderLayout());
        c2.add(eigenLabel, BorderLayout.CENTER);
        c2.add(eigenSlider, BorderLayout.SOUTH);

        Container content = frame.getContentPane();
        content.add(BorderLayout.WEST, c);
        content.add(BorderLayout.EAST, c2);

        EventQueue.invokeLater(()->{
            frame.pack();
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            frame.setVisible(true);
        });

        List<BufferedImage> images = values.stream().map(this::createMap).collect(Collectors.toList());
        slider.addChangeListener(evt->{
            int i = slider.getValue();
            label.setIcon(new ImageIcon(images.get(i)));
            label.setText(i + "/" + slider.getMaximum());
        });

        slider.setValue(0);

        System.out.println("training");
        Trainer train = new Trainer(values);
        train.calculateEigenVectors();
        List<BufferedImage> eigenImages = train.getEigenVectors().stream().map(this::createMap).collect(Collectors.toList());


        eigenSlider.addChangeListener(evt->{
            int i = eigenSlider.getValue();
            eigenLabel.setIcon(new ImageIcon(eigenImages.get(i)));
            eigenLabel.setText(i + "/" + eigenSlider.getMaximum());
        });
        eigenSlider.setValue(0);
    }

    public BufferedImage createMap(double[] values){
        return createMap(width, height, xsize, ysize, values);
    }

    static public BufferedImage createMap(int width, int height, int xsize, int ysize, double[] values){
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;

        for(int i = 0; i<values.length; i++){
            double d = values[i];
            min = d<min?d:min;
            max = d>max?d:max;
        }



        BufferedImage img = new BufferedImage(width*xsize, height*ysize, BufferedImage.TYPE_INT_ARGB);
        Graphics g = img.getGraphics();
        for(int i = 0; i<width; i++){
            for(int j = 0; j<height; j++){

                double v = values[i*width + j]; //only display the first channel.

                if(v<0){
                    Color c = new Color(0, 0, (int)(255*v/min));
                    g.setColor(c);
                } else{
                    Color c = new Color((int)(255*v/max), 0, 0);
                    g.setColor(c);
                }

                //g.fillRect(i*xsize, j*ysize, xsize, ysize);
                g.fillPolygon(
                        new int[]{i*xsize, (i+1)*xsize, i*xsize},
                        new int[]{j*ysize, j*ysize, (j+1)*ysize},
                        3
                );
                /*
                v = values[(i*width + j)*channels + 1]; //only display the first channel.

                if(v<0){
                    Color c = new Color(0, 0, (int)(255*v/min));
                    g.setColor(c);
                } else{
                    Color c = new Color((int)(255*v/max), 0, 0);
                    g.setColor(c);
                }*/
                g.fillPolygon(
                        new int[]{(i+1)*xsize, (i+1)*xsize, i*xsize},
                        new int[]{j*ysize, (j+1)*ysize, (j+1)*ysize},
                        3
                );

            }
        }
        g.dispose();

        return img;
    }

    public static void main(String[] args){

        TwoDHeatMap map = new TwoDHeatMap(32, 32);
        FileDialog log = new FileDialog((JFrame)null, "Choose vector data");
        log.setMode(FileDialog.LOAD);
        log.setVisible(true);
        String f = log.getFile();
        if(f==null) return;
        String d = log.getDirectory();
        try {
            map.displayData(SimpleReader.loadData(Paths.get(d,f)));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
