
#include <iostream>
#include <cmath>

class Field {
protected:
    double *value;

public:
    // constructors and destructors
    Field(){
        value=new double[3]{0.0, 0.0, 0.0};
    }

    Field(double x,double y,double z) {
        value=new double[3]{x,y,z};
    }

    Field(const Field &other){
        value=new double[3]{other.value[0],other.value[1],other.value[2]};
    }

    virtual 
    ~Field(){
        delete[] value;
    }

    virtual 
    void print_mag()const{
        std::cout<<"Components:("<<value[0]<<","<<value[1]<<","<<value[2]<<")"<<std::endl;
    }

    friend 
    std::ostream &operator<<(std::ostream &os, const Field &f){
        os<<"("<<f.value[0]<<","<<f.value[1]<<","<<f.value[2]<<")";
        return os;
    }
};

class ElectricField: 
public Field{
private:
    double calc_electric_field;

public:
    // constructors
    ElectricField(double x=0.0,double y=0.0, double z=0.0) 
        :Field(x,y,z),calc_electric_field(0.0){}

    ElectricField(const ElectricField &other) 
        :Field(other), calc_electric_field(other.calc_electric_field){}

    void calculate_electric_field(double Q,double r){
        const double ep0 = 8.85e-12;
        calc_electric_field=Q/(4*M_PI*r*r*ep0);
    }

    ElectricField operator+(const ElectricField &other)const{
        return ElectricField(value[0]+other.value[0],value[1]+other.value[1],value[2]+other.value[2]);
    }

    friend 
    std::ostream &operator<<(std::ostream &os,const ElectricField &e){
        os <<"Electric Field ";
        e.print_mag();
        os <<"Calculated Electric Field:"<<e.calc_electric_field;
        return os;
    }
};

class MagneticField 
:public Field{
private:
    double calc_magnetic_field;

public:
    //constructors
    MagneticField(double x=0.0,double y=0.0, double z=0.0) 
        :Field(x,y,z),calc_magnetic_field(0.0){}

    MagneticField(const MagneticField &other):Field(other),calc_magnetic_field(other.calc_magnetic_field) {}

    void calculate_magnetic_field(double I,double r){
        const double mu_0=4*M_PI*1e-7;
        calc_magnetic_field=(mu_0*I)/(2*M_PI*r);
    }

    MagneticField operator+(const MagneticField &other) const{
        return MagneticField(value[0]+other.value[0],value[1]+other.value[1],value[2]+other.value[2]);
    }

    friend 
    std::ostream &operator<<(std::ostream &os,const MagneticField &m){
        os <<"Magnetic Field ";
        m.print_mag();
        os << "Calculated Magnetic Field:"<< m.calc_magnetic_field;
        return os;
    }
};

int main() {

    ElectricField eField(0, 1e5, 1e3);
    MagneticField mField(1e-3, 0, 1e-4);

    std::cout<<"Initial Electric Field:\n";
    eField.print_mag();

    std::cout<<"\nInitial Magnetic Field:\n";
    mField.print_mag();

    eField.calculate_electric_field(1e-6, 0.1);  
    mField.calculate_magnetic_field(10, 0.05); 

    std::cout<<"\nCalculated Electric Field:\n"<<eField<<std::endl;
    std::cout<<"\nCalculated Magnetic Field:\n"<<mField<<std::endl;

    ElectricField eField2(1e4,2e4,3e4);
    ElectricField eFieldSum=eField+eField2;
    std::cout<<"\nSum of Electric Fields:\n"<<eFieldSum<<std::endl;

    MagneticField mField2(1e-3, 2e-3, 3e-3);
    MagneticField mFieldSum=mField+mField2;
    std::cout<<"\nSum of Magnetic Fields:\n"<<mFieldSum<<std::endl;

    return 0;
}

