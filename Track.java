package runsimulation1;
// this import is needed for the file input/output functionality
import java.io.*;

class Track
{
    // store momentum components and positions of the Particle
    // mass is a fixed value
    
    float mass;

    float [] px;
    float [] py;
    float [] pz;

    float [] t;
    float [] x;
    float [] y;
    float [] z;

    int lastpos;

    Track(double m, int steps) {
        mass = (float) m;
        px = new float[steps];
        py = new float[steps];
        pz = new float[steps];
        t = new float[steps];
        x = new float[steps];
        y = new float[steps];
        z = new float[steps];
        lastpos = 0;
    }

    void savePositionMomentum(Particle p) {
        if (lastpos >= px.length) {
            System.err.println("Out of space in output track storage, allow for more steps.");
            return;
        }
        px[lastpos] = (float) p.px;
        py[lastpos] = (float) p.py;
        pz[lastpos] = (float) p.pz;
        x[lastpos] = (float) p.x;
        y[lastpos] = (float) p.y;
        z[lastpos] = (float) p.z;
        t[lastpos] = (float) p.t;
        lastpos++;
    }
        
    public void writeToDisk(String filename)
    {
        // this sends the output to a file with name "filename"
        // the block with try { ... } catch (IOException e) { ... } is needed to handle the case,
        // where opening the file fails, e.g. disk is full or similar
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file " + filename + ". Track data was not saved.");
            return;
        }

        // now make a loop to write the contents of the arrays to disk, one number at a time
        outputFile.println("time [s], x [m], y [m], z[m], E [MeV], px [MeV], py [MeV], pz [MeV]");
        for (int n = 0; n < lastpos; n++) {
            // comma separated values
            double E = Math.sqrt(mass*mass + px[n]*px[n] + py[n]*py[n] + pz[n]*pz[n]);
            outputFile.println(t[n] + "," + x[n] + "," + y[n] + "," + z[n] + ","
                               + E + "," + px[n] + "," + py[n] + "," + pz[n]);
        }
        outputFile.close(); // close the output file
        return;
    }
    
}