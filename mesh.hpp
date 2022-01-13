#ifndef mesh_hpp
#define mesh_hpp

namespace dauphine
{

class mesh 
{
    public:
        mesh(spot,dx,T,dt);
    private:
        double n_dt;
        double n_dx;
}

}


#endif