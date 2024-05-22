
struct profile {
   float DeltaDiff;
   float DeltaCum;
   float Velocity;
   float Ri,Rm,Rs;
};

void compute_profiles();
void bin2ascii_profile(int voidID); 
