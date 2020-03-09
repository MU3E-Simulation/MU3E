package runsimulation1;

class Particle
{
    // store the momentum components of the Particle
    // note we store mass rather than energy
    
    double m;
    double px;
    double py;
    double pz;

    // integer charge in units of elementary charge
    int Q;

    // store the time-position components of the Particle
    double t;
    double x;
    double y;
    double z;

    Particle()
    {
        // default constructor: set everything to zero
        m = px = py = pz = 0.;
        t = x = y = z = 0.;
        Q = 0;
    }
    
    Particle(Particle in)
    {
        // constructor that initialises from an existing object
        m = in.m;
        px = in.px;
        py = in.py;
        pz = in.pz;
        t = in.t;
        x = in.x;
        y = in.y;
        z = in.z;
        Q = in.Q;
    }
    
    double momentum()
    {
        double p = px*px + py*py + pz*pz;
        p = Math.sqrt(p);
        return p;
    }

    double mass()
    {
        return m;
    }

    double E()
    {
        double p = momentum();
        return Math.sqrt(m*m + p*p);
    }

    double gamma()
    {
        return (E()/mass());
    }

    double beta()
    {
        return (momentum()/E());
    }

    double radius(double B)
    {
        return (momentum()/(300.*B*Math.abs(Q)));
    }

    void print() {
        System.out.println("Particle with mass = " + m + " MeV, charge = " + Q
                           + ", momentum = " + momentum() + " MeV , beta = " + beta());
        System.out.println("(E, px, py, pz) = (" + E() + ", " + px + ", " + py + ", " + pz + ") MeV");
        System.out.println("(t, x, y, z) = (" + t + ", " + x + ", " + y + ", " + z + ") m");
    }

    public void applySmallRotation(double dtheta_xz, double dtheta_yz)
    {
        // calculate new momentum direction
        // approximates this into two sequential x-z and y-z changes
        double p, theta;
        
        p = Math.sqrt(px*px + pz*pz);
        theta = Math.atan2(px, pz) + dtheta_xz;
        px = p*Math.sin(theta);
        pz = p*Math.cos(theta);

        p = Math.sqrt(py*py + pz*pz);
        theta = Math.atan2(py, pz) + dtheta_yz;
        py = p*Math.sin(theta);
        pz = p*Math.cos(theta);
    }
    
    public void reduceEnergy(double Eloss)
    {
        // reduce energy of particle while keeping direction the same

        double Enew = E() - Eloss;
        if (Enew <= m) {
            // all kinetic energy lost, put particle at rest
            px = 0.;
            py = 0.;
            pz = 0.;
        } else {
            double pnew = Math.sqrt(Enew*Enew - m*m);
            double factor = pnew/momentum();
            px = px * factor;
            py = py * factor;
            pz = pz * factor;
        }
    }
    

    public double [] getState()
    {
        double [] state = {t, x, y, z, E(), px, py, pz};
        return state;
    }

    public void setState(double [] in)
    {
        t = in[0];
        x = in[1];
        y = in[2];
        z = in[3];
        // Energy is NOT set, we assume this is given by sqrt(p^2 + m^2) with m=const.
        px = in[5];
        py = in[6];
        pz = in[7];
    }

    public void setState(Particle p)
    {
        t = p.t;
        x = p.x;
        y = p.y;
        z = p.z;
        m = p.m;
        px = p.px;
        py = p.py;
        pz = p.pz;
    }

    public double distance(Particle p)
    {
        double dx = x-p.x;
        double dy = y-p.y;
        double dz = z-p.z;
        return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }
}