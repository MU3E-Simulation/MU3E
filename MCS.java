package runsimulation1;


class MCS
{
    double rho;
    double Z;
    double A;
    
    
    public MCS(double rho_in, double Z_in, double A_in)
    {
        rho=rho_in;
        Z=Z_in;
        A=A_in;
    }

    public double getX0()
    {
        if (Z==0 ){
            return 0.;
        }
        if (A==0){
            return 0.;
        }
        if (rho==0){
            return 0.;
        }
        double X0 = (716.4*A)/(rho*Z*(Z+1)*Math.log(287/Math.sqrt(Z)));
        // shall return X0 in m
        return X0*0.01;
        
    }

    public double getTheta0(Particle part, double x)
    {
        if (Z==0 ){
            return 0.;
        }
        if (A==0){
            return 0.;
        }
        if (rho==0){
            return 0.;
        }
        double charge= -1;
        double X0 = getX0();
        double beta = part.beta();
        double momentum = part.momentum();
        double theta0 = (13.6/(beta*momentum))*charge*Math.sqrt(x/X0)*(1+0.038*Math.log(x/X0));
               
        // shall return Theta0 for material thickness x
        return theta0;
    }
}