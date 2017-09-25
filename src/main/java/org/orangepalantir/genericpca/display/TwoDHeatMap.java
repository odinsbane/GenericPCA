package org.orangepalantir.genericpca.display;

import org.orangepalantir.genericpca.SimpleReader;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JProgressBar;
import javax.swing.JSlider;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.EventQueue;
import java.awt.FileDialog;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by msmith on 03.08.17.
 */
public class TwoDHeatMap {
    int width;
    int height;
    int channels;
    int xsize = 10;
    int ysize = 10;
    int dex = 0;
    public TwoDHeatMap(int width, int height, int channels){
        this.width = width;
        this.height = height;
        this.channels = channels;
    }
    public void displayData(List<double[]> values){

        JFrame frame = new JFrame("components");

        JLabel label = new JLabel("waiting");
        JSlider slider = new JSlider();
        slider.setMaximum(values.size()-1);
        slider.setMinimum(0);

        Container c = frame.getContentPane();
        c.add(label, BorderLayout.CENTER);
        c.add(slider, BorderLayout.SOUTH);

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

    }
    public BufferedImage createMap(double[] values){
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

                double v = values[(i*width + j)*channels + dex]; //only display the first channel.

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
        TwoDHeatMap map = new TwoDHeatMap(10, 10, 2);
        FileDialog log = new FileDialog((JFrame)null, "Choose vector data");
        log.setMode(FileDialog.LOAD);
        log.setVisible(true);
        String f = log.getFile();
        String d = log.getDirectory();
        try {

            map.displayData(SimpleReader.loadData(Paths.get(d,f)));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
