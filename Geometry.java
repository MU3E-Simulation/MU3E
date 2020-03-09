package runsimulation1;

import java.util.Random;

class Geometry
{
    // Class to define the Geometry of an experiment
    // the "type" is a number defining the basic geometric object
    // the "rho_Z_A" array stores material information: 
    //    density, atomic number (nucleus charge), mass number (protons+neutrons)
    // the "shape" array stores the geometrical information
    // * cuboid: type=1,
    //   shape = 6 values representing two corners of the internal diagonal (x0,y0,z0), (x1,y1,z1)
    // every basic object can be identified by a number
    // the first object is taken as to be the "world",
    //    i.e. we can use it to abort the simulation
    // maxShapes is the maximum allowed number of basic objects, extend if needed

    // the class also impelments thee additional helper functions:
    // * doEloss -- calculate EnergyLoss via associated class and apply to Particle
    // * doMultScatter -- calculate MultipleScattering via associated class and apply to Particle
    // * detectParticle -- simulate a detected particle position
    
    private Random randGen;
    private static final int maxShapes = 100;
    private int nshapes;
    private int [] type;
    private double [] rho;
    private double [] Z;
    private double [] A;

    private double [][] shapes;
    
    private EnergyLoss [] Eloss;
    private MCS [] MultScatter;

    private double minfeaturesize;

    public Geometry(double featuresize)
    {
        randGen = new Random();
        nshapes = 0;
        type = new int[maxShapes];
        rho = new double[maxShapes];
        Z = new double[maxShapes];
        A = new double[maxShapes];
        shapes = new double[maxShapes][];
        
        Eloss = new EnergyLoss[maxShapes];
        MultScatter = new MCS[maxShapes];

        minfeaturesize = featuresize;
    }

    public int getNshapes() { return nshapes; }

    public int AddCuboid(double x0, double y0, double z0,
                         double x1, double y1, double z1,
                         double rhoin, double Zin, double Ain)
    {
        if (nshapes >= maxShapes) {
            return -1;
        }
        
        type[nshapes] = 1;
        shapes[nshapes] = new double[6];
        shapes[nshapes][0] = x0;
        shapes[nshapes][1] = y0;
        shapes[nshapes][2] = z0;
        shapes[nshapes][3] = x1;
        shapes[nshapes][4] = y1;
        shapes[nshapes][5] = z1;

        rho[nshapes] = rhoin;
        Z[nshapes] = Zin;
        A[nshapes] = Ain;

        
        Eloss[nshapes] = new EnergyLoss(rhoin, Zin, Ain);
        MultScatter[nshapes] = new MCS(rhoin, Zin, Ain);

        nshapes++;
        return (nshapes-1);
    }
    
    public void Print()
    {
        System.out.println("stored " + getNshapes() + " objects.");
        for (int i = 0; i < nshapes; i++) {
            if (i == 0) {
                System.out.println("Maximum size of experiment given by object 0:");
            }
            if (type[i] == 1) {
                System.out.println("Geometry object #" + i + " = cuboid.");
                System.out.printf("   corners (%f, %f, %f) - (%f, %f, %f)%n",
                                  shapes[i][0], shapes[i][1], shapes[i][2],
                                  shapes[i][3], shapes[i][4], shapes[i][5]);
            }
            System.out.printf("   material rho = %.3f g/cm^3, Z = %f, A = %f%n",
                              rho[i], Z[i], A[i]);
        }
        System.out.println("When scanning for volume transitions, the smallest feature size discovered will be " + minfeaturesize + "m.");
    }

    public boolean isInVolume(double x, double y, double z, int id)
    {
        // test if point (x,y,z) is in volume with identifier id
        
        // abort if being asked for an undefined volume
        if (id >= getNshapes()) {
            return false;
        }
        
        if (type[id] == 1) {
            // cuboid
            return ( shapes[id][0] <= x
                     && shapes[id][1] <= y
                     && shapes[id][2] <= z
                     && x <= shapes[id][3]
                     && y <= shapes[id][4]
                     && z <= shapes[id][5] );
        }
        
        return false;
    }
    
    public int getVolume(double x, double y, double z)
    {
        // cycle through volumes in opposite order
        for (int i = getNshapes()-1; i >= 0; i--) {
            if (isInVolume(x, y, z, i)) {
                return i;
            }
        }
        
        // if we arrived here, we are outside everything
        return -1;
    }

    public boolean isInVolume(Particle p, int id)
    {
        // test if particle p is currently in volume with identifier id
        return isInVolume(p.x, p.y, p.z, id);
    }
    
    public int getVolume(Particle p)
    {
        // get the highest volume number the particle is in
        return getVolume(p.x, p.x, p.z);
    }
    
    public void doEloss(Particle p, double dist)
    {
        int volume = getVolume(p);
        
        if (volume >= 0) {
            double lostE = Eloss[volume].getEnergyLoss(p)*dist;
            p.reduceEnergy(lostE);
        }
    }
    
    public void doMultScatter(Particle p, double dist)
    {
        int volume = getVolume(p);
        
        if (volume < 0) {
            return;
        }
        
        double theta0 = MultScatter[volume].getTheta0(p, dist);

        if (Math.abs(theta0) > 0.) {
            p.applySmallRotation(randGen.nextGaussian()*theta0,
                                 randGen.nextGaussian()*theta0);
        }
        
    }

    public double [][] detectParticles(Track simtrack)
    {
        double [][] detection = new double[getNshapes()][4];

        // loop over each volume and average over the matching points
        for (int idet = 1; idet < getNshapes(); idet++) {

            // count points in volume
            int ncross = 0;

            for (int ipoint = 0; ipoint < simtrack.lastpos; ipoint++) {
                if (getVolume(simtrack.x[ipoint], simtrack.y[ipoint], simtrack.z[ipoint]) == idet) {
                    detection[idet][0] += simtrack.t[ipoint];
                    detection[idet][1] += simtrack.x[ipoint];
                    detection[idet][2] += simtrack.y[ipoint];
                    detection[idet][3] += simtrack.z[ipoint];
                    ncross++;
                }
            }
            if (ncross > 0) {
                for (int i = 0; i < 4; i++) {
                    detection[idet][i] /= ncross;
                }
            }
        }
        
        return detection;
    }

    public double scanVolumeChange(Particle p, Particle pold, int lastVolume)
    {
        // scan the straight line between last and new position, if we missed a
        // feature of the experiment
        double dist = p.distance(pold);
        int nsteps = (int) (dist/minfeaturesize) + 1;

        if (nsteps <= 1) {
            // step was small enough, continue
            return 1.;
        }
        
        double [] pos = {pold.x, pold.y, pold.z};
        double [] end = {p.x, p.y, p.z};
        double [] delta = new double[3];
        for (int i = 0; i < 3; i++) {
            delta[i] = (end[i]-pos[i])/nsteps;
        }
        for (int n = 1; n <= nsteps; n++) {
            if (getVolume(pos[0]+n*delta[0], pos[1]+n*delta[1], pos[2]+n*delta[2]) != lastVolume) {
                if (n == 1) {
                    return (0.5 / nsteps);
                } else {
                    return (n-1.0) / nsteps;
                }
            }
        }
        return 1.;
    }
}