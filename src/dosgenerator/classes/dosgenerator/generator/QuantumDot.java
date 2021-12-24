/*
 * Copyright (C) 2020-2021 Alban Lafuente
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
import com.github.kilianB.pcg.fast.PcgRSFast;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.nevec.rjm.BigDecimalMath;

/**
 *
 * @author Alban Lafuente
 */
public class QuantumDot
{
    private final ArrayList<BigDecimal> m_listOfStates = new ArrayList<>();
    private final BigDecimal m_positionX;
    private final BigDecimal m_positionY;
    private final BigDecimal m_radius;
    private final BigDecimal m_height;
    private final int m_numberOfStates;
    private final HashMap<Double, BigDecimal> m_probabilitiesPerlevel;
    private final TreeSet<Double> m_recombinationProbaTree;
    
    private int m_numberOfFreeStates;
    
    public QuantumDot (BigDecimal p_positionX, BigDecimal p_positionY, BigDecimal p_radius, BigDecimal p_height, Map<Double, BigDecimal> p_energyLevelsPopProba, Set<Double> p_recombProba, int p_nbLevels, int p_nbFreeLevels)
    {
        m_positionX = new BigDecimal(p_positionX.toString());
        m_positionY = new BigDecimal(p_positionY.toString());
        m_radius = new BigDecimal(p_radius.toString());
        m_height = new BigDecimal(p_height.toString());
        m_numberOfStates = p_nbLevels;
        m_numberOfFreeStates = p_nbFreeLevels;
        
        m_probabilitiesPerlevel = new HashMap<>();
        for (Double proba: p_energyLevelsPopProba.keySet())
        {
            m_probabilitiesPerlevel.put(proba, new BigDecimal(p_energyLevelsPopProba.get(proba).toString()));
        }
        
        m_recombinationProbaTree = new TreeSet<>();
        for (Double proba: p_recombProba)
        {
            m_recombinationProbaTree.add(proba);
        }
    }

    //ΔEg(InAs/GaAs) ~ 1.1 eV
    public QuantumDot (BigDecimal p_positionX, BigDecimal p_positionY, BigDecimal p_radius, BigDecimal p_height, Metamaterial p_sampleMaterial)
    {
        BigDecimal two = new BigDecimal("2");
        Material QDMaterial = p_sampleMaterial.getMaterial("QD");
        Material barrierMaterial = p_sampleMaterial.getMaterial("barrier");
        
        m_positionX = p_positionX;
        m_positionY = p_positionY;
        
        m_radius = p_radius.multiply(BigDecimal.ONE);
        m_height = p_height.multiply(BigDecimal.ONE);
        BigDecimal equivalentSquareSide = m_radius.multiply(BigDecimalMath.sqrt(BigDecimalMath.pi(MathContext.DECIMAL128), MathContext.DECIMAL128));
        
        BigDecimal CBOffset = p_sampleMaterial.getOffset(QDMaterial.getMaterialName(), barrierMaterial.getMaterialName()); //from https://aip.scitation.org/doi/abs/10.1063/1.125965

        //calculating hole confinement energy, only considering one level
        BigDecimal VBOffset = barrierMaterial.getBandgap().subtract(QDMaterial.getBandgap()).subtract(CBOffset);
        BigDecimal planeEnergyParameterHole = energyParameter(0, equivalentSquareSide, VBOffset, QDMaterial.getHoleEffectiveMass());
        BigDecimal heightEnergyParameterHole = energyParameter(0, m_height, VBOffset, QDMaterial.getHoleEffectiveMass());
        BigDecimal holeConfinementEnergy = (two.multiply(PhysicsVariables.hbar.pow(2)).divide(QDMaterial.getHoleEffectiveMass(), MathContext.DECIMAL128)).multiply(heightEnergyParameterHole.add(two.multiply(planeEnergyParameterHole)));
        
        int nbStates = 0;
        TreeSet<BigDecimal> energyLevels = new TreeSet<>();
        for (int nz = 0 ; nz < 10 ; nz += 1)
        {
            for (int nx = 0 ; nx < 100 ; nx += 1)
            {
                for (int ny = 0 ; ny < 100 ; ny += 1)
                {
                    BigDecimal xEnergyParameterElectron = energyParameter(nx, equivalentSquareSide, CBOffset, QDMaterial.getElectronEffectiveMass());
                    BigDecimal yEnergyParameterElectron = energyParameter(ny, equivalentSquareSide, CBOffset, QDMaterial.getElectronEffectiveMass());
                    BigDecimal zEnergyParameterElectron = energyParameter(nz, p_height, CBOffset, QDMaterial.getElectronEffectiveMass());
                    
                    if (xEnergyParameterElectron.compareTo(BigDecimal.ZERO) < 0 || yEnergyParameterElectron.compareTo(BigDecimal.ZERO) < 0 || zEnergyParameterElectron.compareTo(BigDecimal.ZERO) < 0)
                    {
                        break;
                    }
                    
                    BigDecimal energyXElectron = xEnergyParameterElectron.divide(equivalentSquareSide, MathContext.DECIMAL128).pow(2);
                    BigDecimal energyYElectron = yEnergyParameterElectron.divide(equivalentSquareSide, MathContext.DECIMAL128).pow(2);
                    BigDecimal energyZElectron = zEnergyParameterElectron.divide(p_height, MathContext.DECIMAL128).pow(2);

                    BigDecimal electronConfinementEnergy = (two.multiply(PhysicsVariables.hbar.pow(2)).divide(QDMaterial.getElectronEffectiveMass(), MathContext.DECIMAL128)).multiply(energyXElectron.add(energyYElectron).add(energyZElectron));
                    if (electronConfinementEnergy.compareTo(CBOffset) > 0)
                    {
                        break;
                    }
                    
                    energyLevels.add(electronConfinementEnergy);
                    nbStates += 2;
                    
                    //adding the energy to the list of states
                    BigDecimal totalRecombinationEnergy = electronConfinementEnergy.add(QDMaterial.getBandgap()).add(holeConfinementEnergy);
                    if (totalRecombinationEnergy.compareTo(BigDecimal.ZERO) < 0)
                    {
                        throw new InternalError("Negative recombination energy.");
                    }
                    m_listOfStates.add(totalRecombinationEnergy);
                    m_listOfStates.add(totalRecombinationEnergy);
                }
            }
        }
        
        m_numberOfStates = nbStates;
        m_numberOfFreeStates = m_numberOfStates;
        
        m_probabilitiesPerlevel = new HashMap<>();
        m_recombinationProbaTree = new TreeSet<>();
        
        if (nbStates != 0)
        {
            /**RECOMB PROBA PER LEVEL
             * calculate probability for each level using Fermi-Dirac distribution and the energy calculated from the QD material CB position
             * BIG approximation: chemical potential = 0
             */

            //calcul of the probability for an electron to be on a given level when it recombines
            BigDecimal sumOfProba = BigDecimal.ZERO;
            HashMap<BigDecimal, BigDecimal> levelsProbabilities = new HashMap<>();
            for (BigDecimal energy: energyLevels)
            {
                BigDecimal fermiDiracProba = BigDecimal.ONE.divide(BigDecimal.ONE.add(BigDecimalMath.exp((energy.subtract(BigDecimal.ZERO)).divide(PhysicsVariables.KB.multiply(new BigDecimal("300")), MathContext.DECIMAL128), MathContext.DECIMAL128)), MathContext.DECIMAL128);
                sumOfProba = sumOfProba.add(fermiDiracProba);
                levelsProbabilities.put(energy, fermiDiracProba);
            }

            //normalizing the probabilities, taking care of it reaching all the way to 1, and saving the energy as the complete recombination energy (considering recombination occur toward the highest hole energy) and not the confinement energy 
            BigDecimal sumOfPreviousProba = BigDecimal.ZERO;
            BigDecimal meanEnergy = BigDecimal.ZERO;
            for (BigDecimal energy: energyLevels)
            {
                BigDecimal normalizedProba = levelsProbabilities.get(energy).divide(sumOfProba, MathContext.DECIMAL128);
                sumOfPreviousProba = sumOfPreviousProba.add(normalizedProba);
                if (energy.compareTo(energyLevels.last()) == 0  && sumOfPreviousProba.compareTo(BigDecimal.ONE) != 0)
                {
                    sumOfPreviousProba = BigDecimal.ONE;
                }

                BigDecimal totalRecombinationEnergy = energy.add(QDMaterial.getBandgap()).add(holeConfinementEnergy);
                meanEnergy = meanEnergy.add(totalRecombinationEnergy.multiply(normalizedProba));
                m_recombinationProbaTree.add(sumOfPreviousProba.doubleValue());
                m_probabilitiesPerlevel.put(sumOfPreviousProba.doubleValue(), totalRecombinationEnergy);
            }
        }
    }
    
    public QuantumDot copy()
    {
        return new QuantumDot(m_positionX, m_positionY, m_radius, m_height, m_probabilitiesPerlevel, m_recombinationProbaTree, m_numberOfStates, m_numberOfFreeStates);
    }
    
    public QuantumDot copyWithSizeChange(BigDecimal p_sizeMultiplier, Metamaterial p_sampleMaterial)
    {
        BigDecimal newRadius = m_radius;
        BigDecimal newHeight = m_height;
        
        if (m_radius.compareTo(m_height) > 0)
        {
            newHeight = newHeight.multiply(p_sizeMultiplier);
        }
        else
        {
            newRadius = newRadius.multiply(p_sizeMultiplier);
        }
        
        return new QuantumDot(m_positionX, m_positionY, newRadius, newHeight, p_sampleMaterial);
    }
    
    /**
     * See https://en.wikipedia.org/wiki/Finite_potential_well
     * @param index 
     * @param size
     * @param bandOffset 
     * @param effectiveMass 
     * @return 
     */
    private BigDecimal energyParameter (int index, BigDecimal size, BigDecimal bandOffset, BigDecimal effectiveMass)
    {
        double u02 = (effectiveMass.multiply(size.pow(2)).multiply(bandOffset).divide((new BigDecimal(2)).multiply(PhysicsVariables.hbar.pow(2)), MathContext.DECIMAL128)).doubleValue();
        double vi = 0;
        
        //vi has to be between i*pi/2 and (i+1)*v/2. Minimum Vi should also always be lower than u0
        double minVi = index * Math.PI/2;
        if (Math.pow(minVi, 2) >= u02)
        {
//            System.err.println("Too high of a minVi");
            vi = -1;
        }
        else
        {
            vi = minVi + Math.random()*Double.min(Math.PI/2, Math.sqrt(u02) - minVi);
            
            double maxVi = (index + 1) * Math.PI/2;
            double error = 1E-14;
            double epsilon = 1E-15;
            int counter = 0;
            
            do
            {
                double derivative = derivativeFunction(index, vi);
                if (Math.abs(derivative) <= epsilon)
                {
                    break;
                }
                
                vi = Math.abs(vi - ((functionToOptimize(index, vi) - u02) / derivative));
                
                while (vi <= minVi || vi >= maxVi)
                {
                    //vi has to be between i*pi/2 and (i+1)*pi/2
                    if (vi < minVi)
                    {
                        vi = vi - (Math.PI/2) * (int) ((vi)/(Math.PI/2)) + minVi ;
                    }
                    else
                    {
                        if (vi > maxVi)
                        {
                            vi = vi - (Math.PI/2) * (int) ((vi)/(Math.PI/2)) + minVi;
                        }
                        else
                        {
                            vi *= 1.1;
                        }
                    }
                }
                
                counter += 1;
                if (counter%100 == 0)
                {
                    error *= 2;
                }
            }while(Math.abs(functionToOptimize(index, vi) - u02) >= error);
        }
        
        return new BigDecimal(vi);
    }
    
    private double functionToOptimize(int index, double v)
    {
        if (index % 2 == 0)
        {
            return Math.pow(v, 2) * (1 + Math.pow(Math.tan(v), 2));
        }
        else
        {
            return Math.pow(v, 2) * (1 + 1 / Math.pow(Math.tan(v), 2));
        }
    }
    
    private double derivativeFunction(int index, double v)
    {
        double function = 0;
        double modif = 0;
        
        if (index % 2 == 0)
        {
            function = Math.tan(v);
            modif = 1 / Math.pow(Math.cos(v), 2);
        }
        else
        {
            function = 1 / Math.tan(v);
            modif = 1 / Math.pow(Math.sin(v), 2);
        }
        
        return 2 * v * (1 + Math.pow(function, 2) + v * function * modif);
    }
    
    public BigDecimal getRadius()
    {
        return m_radius;
    }
    
    public ArrayList<BigDecimal> getStates()
    {
        ArrayList<BigDecimal> listOfStates = new ArrayList<>();
        
        for (BigDecimal state: m_listOfStates)
        {
            listOfStates.add(new BigDecimal(state.toPlainString()));
        }
        
        return listOfStates;
    }
    
    public String scaledString(BigDecimal p_sizeScale)
    {
        //new scale: number.scale() - number.precision() gives the number of digits after the point in scientific notation. Setting the scale to this + 11 gives us at least 10 digits after the points, which is enough
        BigDecimal scaledX = (m_positionX.divide(p_sizeScale, MathContext.DECIMAL128)).setScale(m_positionX.scale() - m_positionX.precision() + 11, RoundingMode.HALF_UP).stripTrailingZeros();
        BigDecimal scaledY = (m_positionY.divide(p_sizeScale, MathContext.DECIMAL128)).setScale(m_positionY.scale() - m_positionY.precision() + 11, RoundingMode.HALF_UP).stripTrailingZeros();
        BigDecimal scaledRadius = (m_radius.divide(p_sizeScale, MathContext.DECIMAL128)).setScale(m_radius.scale() - m_radius.precision() + 11, RoundingMode.HALF_UP).stripTrailingZeros();
        BigDecimal scaledHeight = (m_height.divide(p_sizeScale, MathContext.DECIMAL128)).setScale(m_height.scale() - m_height.precision() + 11, RoundingMode.HALF_UP).stripTrailingZeros();
//        BigDecimal scaledEnergy = (m_energyLevelPopulatedProbabilities.divide(p_energyScale, MathContext.DECIMAL128)).setScale(m_energyLevelPopulatedProbabilities.scale() - m_energyLevelPopulatedProbabilities.precision() + 11, RoundingMode.HALF_UP).stripTrailingZeros();
        
        return scaledX + "\t" + scaledY + "\t" + scaledRadius + "\t" + scaledHeight/* + "\t" + scaledEnergy*/;
    }
    
    @Override
    public String toString()
    {
        return m_positionX + "\t" + m_positionY + "\t" + m_radius + "\t" + m_height + "\t" + m_probabilitiesPerlevel;
    }
    
    public BigDecimal getDistance (BigDecimal p_positionX, BigDecimal p_positionY)
    {
        BigDecimal squaredDistance = ((m_positionX.subtract(p_positionX)).pow(2)).add(((m_positionY.subtract(p_positionY)).pow(2)));
        BigDecimal distance;
        
        if (squaredDistance.compareTo(BigDecimal.ZERO) == 0)
        {
            distance = BigDecimal.ZERO;
        }
        else
        {
            distance = BigDecimalMath.sqrt(squaredDistance);
        }
        
        return distance;
    }
    
    public BigDecimal getX()
    {
        return m_positionX;
    }
    
    public BigDecimal getY()
    {
        return m_positionY;
    }
}
