
struct grid {
  int NumMem;
  int *Member;
  int *Neighbour;
};

struct query {
   vector <int> i;	
   vector <int> j;	
   vector <int> k;
};

template <class T>
void build_grid_list(T Points, int NumPoints, struct grid *GridList, int NumGrid, double *GridSize, bool compute_neigh);
void query_grid(struct query *Query, double *GridSize, double R1, double R2);
void free_query_grid(struct query *Query);
void free_grid_list(struct grid *GridList, int NG);
double max_grid_size(double *GridSize);

