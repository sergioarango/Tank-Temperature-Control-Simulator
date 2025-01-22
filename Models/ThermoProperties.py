from scipy.integrate import quad

class ThermoProperties:
    
    def __init__(self, name = "TermoProperties"):
        return

    def Hliq (self, T, P, Patm):
        self.T=T
        self.P=P
        self.Patm = Patm
        
        self.A=7.243E+01
        self.B=1.039E-02
        self.C=-1.497E-06
        
        self.MWH20 = 18.015  # Water molecular Weight [kg/kmol]
        self.TcH2O = 647.19  #Critcal temperature [K]
        self.PcH2O = 22055   #Critical pressure [kPa]
        self.wH2O = 0.34486  #acentric factor
        self.R = 8.314       #Universal gas constant [J/molK]
        
        #Reference temperature
        self.TrH2O = self.T/self.TcH2O
        
        #Reference Pressure
        self.PrH2O = (self.P+self.Patm)/self.PcH2O
        
        #compressibility factor
        
        #Bo Parameter
        self.BoH2O = 0.083-(0.422/pow(self.TrH2O, 1.6))
        
        #B1 Parameter
        self.B1H2O = 0.139-(0.172/pow(self.TrH2O, 4.2))
        
        #Compressibility factor estimation
        self.ZH2O = 1+(self.BoH2O+self.wH2O*self.B1H2O)*(self.PrH2O/self.TrH2O)
        
        #specific vlolume
        self.vH2O = (self.ZH2O*self.R*self.T)/(self.P+self.Patm)  # [m3/kmol]
        
        def integrate(A, B, C, Tref, T):
                result = quad(lambda x: A + B*x + C*x**2, Tref, T)
                I = (result[0])
                return I
        
        self.Tref = 273.15
        
        MolarEnthalpy = integrate(A=self.A, B=self.B, C=self.C, Tref=self.Tref, T=self.T ) + (self.P+self.Patm)*(self.vH2O/self.MWH20)
        
        return MolarEnthalpy