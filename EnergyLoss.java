package runsimulation1;

class EnergyLoss
{
    double rho;
    double Z;
    double A;
   
    public EnergyLoss(double rho_in, double Z_in, double A_in)
    {
        rho=rho_in;
        Z=Z_in;
        A=A_in;
               
    }
    
    public double getEnergyLoss(Particle p)
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
        double beta = p.beta();
        double gamma = p.gamma();
        double mass = p.mass();
        double K = 0.307075; // MeV cm^2
        double charge = -1;
        double Me = 0.511; // MeV
        
        double wmax = (2 * Me * gamma * gamma * beta * beta) /(1+(2*gamma*Me)/(mass+((Me/mass)*(Me/mass))));
        double I = 0.0000135 * Z;
        double loss = K*charge*charge*rho*(Z/A)*(1/(beta*beta))*(0.5*Math.log(2*Me*beta*beta*gamma*gamma*wmax/(I*I))-(beta*beta));
                
              
        return loss * 100;        
         // shall return energy loss in MeV/m
        
    }
    
}