package runsimulation1;

/// this import is needed for the file input/output functionality
import java.io.*;

class Histogram
{       
    private double binlow, binhigh;
    private double binwidth;
    private int nbins;
    private double[] binCentre;
    private String histname;

    // double array to store the actual histogram data
    private double[] sumWeights;

    private long underflow, overflow;
    private long nfilled;

    // constructor for the class Histogram
    public Histogram(int numberOfBins, double start, double end, String name)
    {
        // store the parameters and setup the histogram
        // note that parameters need to have different names than class variables
        nbins = numberOfBins;
        binlow = start;
        binhigh = end;
        histname = name;

        binwidth = (binhigh - binlow) / (double) nbins;
        sumWeights = new double[nbins];
        underflow = 0;
        overflow = 0;
        nfilled = 0;
        
        // calculate and save the x coordinate of the centre of each bin
        binCentre = new double[nbins];
        for (int i = 0; i < nbins; i++) {
            binCentre[i] = binlow + (i+0.5)*binwidth;
        }
    }

    public int getNbins()
    {
        return nbins;
    }

    public long getUnderflow()
    {
        return underflow;
    }

    public long getOverflow()
    {
        return overflow;
    }
    public long getNfilled()
    {
        return nfilled;
    }

    public void fill(double value)
    {
        if (value < binlow) {
            underflow++;
        } else if (value >= binhigh) {
            overflow++;
        } else {
            // add weight to the correct bin
            int ibin = (int) ( (value - binlow)/binwidth);
            sumWeights[ibin] = sumWeights[ibin] + 1.0;
        }
        nfilled++;
    }

    public double getContent(int nbin)
    {
        // returns the contents on bin 'nbin' to the user
        return sumWeights[nbin];
    }

    public double getError(int nbin)
    {
        // returns the error on bin 'nbin' to the user
        return Math.sqrt(sumWeights[nbin]);
    }
    
    //-------------------------------------
    public void print()
    {
        for (int bin = 0; bin < getNbins(); bin++) {
            System.out.println("Bin " + bin + " = " +getContent(bin)
                               + " +- " + getError(bin));
        }
        System.out.println("The number of fills = " + getNfilled());
        System.out.println("Underflow = " + getUnderflow()
                           + ", Overflow = " + getOverflow());
    }
    
    //-------------------------------------
    public void writeToDisk(String filename)
    {
        // this sends the output to a file with name "filename"
        // the block with try { ... } catch (IOException e) { ... } is needed to handle the case,
        // where opening the file fails, e.g. disk is full or similar
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file " + filename + ". Histogram data was not saved.");
            return;
        }

        // Write the file as a comma seperated file (.csv) so it can be read it into EXCEL
        // first some general information about the histogram
        outputFile.println("histname, " + histname);
        outputFile.println("binlow, " + binlow);
        outputFile.println("binwidth, " + binwidth);
        outputFile.println("nbins, " + nbins);
        outputFile.println("underflow, " + underflow);
        outputFile.println("overflow, " + overflow);

        // now make a loop to write the contents of each bin to disk, one number at a time
        // together with the x-coordinate of the centre of each bin.
        for (int n = 0; n < nbins; n++) {
            // comma separated values
            outputFile.println(n + "," + binCentre[n] + "," + getContent(n) + "," + getError(n));
        }
        outputFile.close(); // close the output file
    }
}