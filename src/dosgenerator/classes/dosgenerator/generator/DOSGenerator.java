/*
 * Copyright (C) 2021 audreyazura
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package dosgenerator.generator;

import albanlafuente.physicstools.physics.Material;
import albanlafuente.physicstools.physics.Metamaterial;
import albanlafuente.physicstools.physics.PhysicsVariables;
import java.math.BigDecimal;
import java.math.MathContext;
import com.github.kilianB.pcg.fast.PcgRSFast;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 *
 * @author audreyazura
 */
public class DOSGenerator
{
    public static void main(String[] args)
    {
        List<QuantumDot> QDList = new ArrayList<>();
        
        /**********************************************************************
         *                      GETTING MATERIALS                             * 
         **********************************************************************/
        
        Map<String, Material> materialMap = new HashMap<>(); 
        Properties InAsProperties = new Properties();
        Properties GaAsProperties = new Properties();
        Properties metamaterialProperties = new Properties();
        SCSVLoader functionLoader = new SCSVLoader();
        try
        {
            InAsProperties.load(new FileReader("ressources/materials/InAs.mat"));
            GaAsProperties.load(new FileReader("ressources/materials/GaAs.mat"));
            metamaterialProperties.load(new FileReader("ressources/metamaterials/InAsGaAs.metamat"));
        }
        catch (IOException ex)
        {
            Logger.getLogger(DOSGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }
        materialMap.put(InAsProperties.getProperty("name"), new Material(InAsProperties, functionLoader));
        materialMap.put(GaAsProperties.getProperty("name"), new Material(GaAsProperties, functionLoader));
        Metamaterial sampleMaterial = new Metamaterial(metamaterialProperties, materialMap);
        
        /**********************************************************************
         *                  LOADING ALREADY GENERATED QDS                     * 
         **********************************************************************/
        
        System.out.println("Loading already generated QDs");
        
        String QDListFile = "/home/audreyazura/Documents/Work/Simulation/DOSEvolution/QDList.dat";
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(QDListFile));
            Pattern numberRegex = Pattern.compile("^\\-?\\d+(\\.\\d+(e(\\+|\\-)\\d+)?)?");

            String line;
            while (((line = fileReader.readLine()) != null))
            {	    
                String[] lineSplit = line.strip().split("[,;\t]");

                if(numberRegex.matcher(lineSplit[0]).matches())
                {
                    BigDecimal x = formatBigDecimal((new BigDecimal(lineSplit[0].strip())).multiply(PhysicsVariables.UnitsPrefix.NANO.getMultiplier()));
                    BigDecimal y = formatBigDecimal((new BigDecimal(lineSplit[1].strip())).multiply(PhysicsVariables.UnitsPrefix.NANO.getMultiplier()));
                    BigDecimal radius = formatBigDecimal(((new BigDecimal(lineSplit[2].strip()))).multiply(PhysicsVariables.UnitsPrefix.NANO.getMultiplier()));
                    BigDecimal height = formatBigDecimal((new BigDecimal(lineSplit[3].strip())).multiply(PhysicsVariables.UnitsPrefix.NANO.getMultiplier()));

                    QuantumDot currentQD = new QuantumDot(x, y, radius, height, sampleMaterial);
                    QDList.add(currentQD);
                }
            }
        }
        catch (IOException ex)
        {
            System.out.println("No QD file found, continuing on full randomized QDs");
        }
        
        /**********************************************************************
         *                     GENERATING THE QD LIST                         * 
         **********************************************************************/
        
        System.out.println("Generating the other needed QDs");
        
        int totalWishedQDs = 100000;
        int nQDs = totalWishedQDs - QDList.size();
        if (nQDs < 0)
        {
            nQDs = 0;
        }
        
        BigDecimal sampleXSize = BigDecimal.ONE.multiply(PhysicsVariables.UnitsPrefix.CENTI.getMultiplier());
        BigDecimal sampleYSize = BigDecimal.ONE.multiply(PhysicsVariables.UnitsPrefix.CENTI.getMultiplier());
        PcgRSFast RNGenerator = new PcgRSFast();
        
        BigDecimal three = new BigDecimal("3");
        for (int i = 0 ; i < nQDs ; i += 1)
        {
            BigDecimal x, y, radiusNano, radius, height;
            QuantumDot createdQD;

            do
            {
                x = formatBigDecimal((new BigDecimal(RNGenerator.nextDouble())).multiply(sampleXSize));
                y = formatBigDecimal((new BigDecimal(RNGenerator.nextDouble())).multiply(sampleYSize));

                do
                {
                    radiusNano = new BigDecimal(RNGenerator.nextGaussian() * 2.1 + 12);
                    radius = formatBigDecimal(radiusNano.multiply(PhysicsVariables.UnitsPrefix.NANO.getMultiplier()));

                }while (radius.compareTo(BigDecimal.ZERO) <= 0);

                do
                {
                    /**
                     * the height is correlated to the radius with the relation
                     * height = radius / 3 - 1.5
                     * with a variation of about +/- 0.5 around the line. To reproduce that variation, we use the relation
                     * height = radius / 3 + GaussianRNG*0.5 - 1.5
                     * GaussianRNG giving a number on a gaussian centered on 0 with a variance of 1.
                     */
                    height = formatBigDecimal((radiusNano.divide(three, MathContext.DECIMAL128)).add(new BigDecimal(RNGenerator.nextGaussian()*0.5 - 1.5)).multiply(PhysicsVariables.UnitsPrefix.NANO.getMultiplier()));
                }while(height.compareTo(BigDecimal.ZERO) <= 0);

                createdQD = new QuantumDot(x, y, radius, height, sampleMaterial);

            }while(!validPosition(createdQD, QDList));

            QDList.add(createdQD);
            System.out.println(QDList.size());
        }
        
        /**********************************************************************
         *                     CALCULATING THE DOS                            * 
         **********************************************************************/
        
        System.out.println("Calculating the DOS");
        
        ArrayList<BigDecimal> everyStates = new ArrayList<>();
        for (QuantumDot QD: QDList)
        {
            everyStates.addAll(QD.getStates());
        }
        
        everyStates.sort(null);
        BigDecimal minState = everyStates.get(0);
        BigDecimal maxState = everyStates.get(everyStates.size() - 1);
        BigDecimal DOSInterval = (new BigDecimal("0.002")).multiply(PhysicsVariables.EV);
        BigDecimal sampleVolume = sampleXSize.multiply(sampleYSize);
        Map<BigDecimal, BigDecimal> DOS = new HashMap<>();
        for (BigDecimal lowestBound = minState ; lowestBound.compareTo(maxState) == -1 ; lowestBound = lowestBound.add(DOSInterval))
        {
            BigDecimal currentMax = lowestBound.add(DOSInterval);
            int nLevels = 0;
            
            while (everyStates.size() > 0 && everyStates.get(0).compareTo(currentMax) <= 0)
            {
                nLevels += 1;
                everyStates.remove(0);
            }
            
            DOS.put(lowestBound, (new BigDecimal(nLevels)).divide(sampleVolume));
        }
        
        /**********************************************************************
         *                       SAVING TO FILES                              * 
         **********************************************************************/
        
        System.out.println("Saving to files");
        
        String DOSDatFile = "/home/audreyazura/Documents/Work/Simulation/DOSEvolution/DOS_" + totalWishedQDs + "QDs.dat";
        
        try
        {
            //writing DOS
            Set<BigDecimal> statesSet = new TreeSet<>(DOS.keySet());
            BufferedWriter DOSwriter = new BufferedWriter(new FileWriter(DOSDatFile));
            DOSwriter.write("Energy (eV)\tDOS (m^-2)");
            for (BigDecimal state: statesSet)
            {
                BigDecimal stateToWrite = state.divide(PhysicsVariables.EV, MathContext.DECIMAL128).setScale(state.scale() - state.precision() + 4, RoundingMode.HALF_UP);

                DOSwriter.newLine();
                DOSwriter.write(stateToWrite.toPlainString() + "\t" + DOS.get(state));
            }
            DOSwriter.flush();
            DOSwriter.close();
            
            BufferedWriter QDWriter = new BufferedWriter(new FileWriter(QDListFile));
            QDWriter.write("x (nm)\ty (nm)\tradius (nm)\theight (nm)");
            for (QuantumDot qd: QDList)
            {
                QDWriter.newLine();
                QDWriter.write(qd.scaledString(PhysicsVariables.UnitsPrefix.NANO.getMultiplier()));
            }
            QDWriter.flush();
            QDWriter.close();
        } 
        catch (IOException ex)
        {
            Logger.getLogger(DOSGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        /**********************************************************************
         *                     MAKING DOS PICTURE                             * 
         **********************************************************************/
        
        String DOSPicFile = "/home/audreyazura/Documents/Work/Simulation/DOSEvolution/DOS_" + totalWishedQDs + "QDs.png";
        
        try
        {
            BufferedWriter gnuplotWriter = new BufferedWriter(new FileWriter("/home/audreyazura/Documents/Work/Simulation/DOSEvolution/.gnuplotScript.gp"));
            gnuplotWriter.write("set terminal png");
            gnuplotWriter.newLine();
            gnuplotWriter.write("set output \"" + DOSPicFile + "\"");
            gnuplotWriter.newLine();
            gnuplotWriter.write("set label \"#QDs: " + totalWishedQDs + "\" at graph 0.02,0.94");
            gnuplotWriter.newLine();
            gnuplotWriter.write("set xlabel \"Energy (eV)\"");
            gnuplotWriter.newLine();
            gnuplotWriter.write("set ylabel \"Density of states (m^{-2})\"");
            gnuplotWriter.newLine();
            gnuplotWriter.write("plot[0.6:1.1] \"" + DOSDatFile + "\" u 1:2 w line notitle");
            gnuplotWriter.newLine();
            gnuplotWriter.write("unset output");
            gnuplotWriter.flush();
            gnuplotWriter.close();

            Runtime.getRuntime().exec("gnuplot /home/audreyazura/Documents/Work/Simulation/DOSEvolution/.gnuplotScript.gp").waitFor();
            Runtime.getRuntime().exec("rm /home/audreyazura/Documents/Work/Simulation/DOSEvolution/.gnuplotScript.gp");
        } 
        catch (IOException|InterruptedException ex)
        {
            Logger.getLogger(DOSGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private static BigDecimal formatBigDecimal(BigDecimal p_toFormat)
    {
        return p_toFormat.stripTrailingZeros();
    }
    
    private static boolean validPosition(QuantumDot p_testedQD, List<QuantumDot> p_existingQDs)
    {
        boolean valid = true;
        
        for (QuantumDot QD: p_existingQDs)
        {
            valid &= p_testedQD.getRadius().add(QD.getRadius()).compareTo(p_testedQD.getDistance(QD.getX(), QD.getY())) < 0;
        }
        
        return valid;
    }
}
