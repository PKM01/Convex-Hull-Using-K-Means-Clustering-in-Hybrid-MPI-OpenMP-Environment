#include <iostream>
#include <cstdlib>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <fstream>

#include "kMeans1.cpp"
#include "con_hull.cpp"

using namespace std;

int main(int argc,char** argv)
{
    double t1,t2;
    int n=15000,k=2; //no of cluster, and total no of point respectively
    int x[n],y[n],convx[n],convy[n],kc[n];
    int rank, size;

    t1=omp_get_wtime();
    MPI_Init(&argc, &argv); //MPI_COMM_WORLD created
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(0));
    if(rank==0)
    {
    // cout<<"Generating Random Numbers...!!!"<<endl;
    // for(int i=0;i<n;i++){
    // 	x[i]=rand()%100;
    // 	y[i]=rand()%100;
    // }
    cout<<"Reading the data points...!!!"<< n << " " << k <<  endl;
    ifstream infile("points_" + to_string(n) + ".txt");
    if (!infile)
    {
        cerr << "Error: Could not open file" << endl;
        return 1;
    }
    for (int i = 0; i < n; i++)
    {
        infile >> x[i] >> y[i];
    }
    }

    MPI_Bcast(x, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, n, MPI_INT, 0, MPI_COMM_WORLD);
    
    //call the clusters' function and get the clusters
    cout<<"Obtaining the clusters...!!"<<endl;
    find_clusters(n,k,x,y,kc);
    MPI_Bcast(kc, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    int yc[n],xc[n],clust_size=0;
    cout<<"rank: "<<rank<<" size: "<<size<<"\n";
    	for(int i=0;i<n;i++){ //put point that belong to same cluster together
       		if(kc[i]==rank){
    			clust_size+=1;
    			xc[clust_size-1]=x[i];
    			yc[clust_size-1]=y[i];
    		}
    	}
    
        MPI_Barrier(MPI_COMM_WORLD);
    
    cout<<"Finding the convex hulls for each cluster obtained..!!"<<endl;
    //for each cluster find the convex hull and store the points
    int con_x[n],con_y[n],con_n=0,*fin_x,*fin_y,fin_n=0;
    //printf("cluster size:\n",clust_size);

    findHull(xc,yc,clust_size,&con_x[0],&con_y[0],&con_n);
    int recv_counts[k];

    MPI_Gather(&con_n, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD); //recv_count will gather no of points in each c_hull to root for final hull calculation 
    MPI_Reduce(&con_n, &fin_n, 1, MPI_INT, MPI_SUM, 0,  MPI_COMM_WORLD); //fin_n contain total no of points for final c_hull calculation
    int recvx[fin_n],recvy[fin_n];

    int displ[k];
        if(rank==0)
        {
            displ[0]=0;
            int sum=recv_counts[0];
            for(int i=1;i<k;i++){
                displ[i]=displ[i-1]+recv_counts[i-1];
                sum+=recv_counts[i];}

        }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gatherv(con_x, con_n, MPI_INT,recvx, recv_counts, displ, MPI_INT,0, MPI_COMM_WORLD);
    MPI_Gatherv(con_y, con_n, MPI_INT,recvy, recv_counts, displ, MPI_INT,0, MPI_COMM_WORLD);


    if(rank==0)
        {
            printf("total:%d\n",fin_n);

            findHull(recvx,recvy,fin_n,&con_x[0],&con_y[0],&con_n);

            printf("final hull:\n");
            for(int i=0;i<con_n;i++)
                printf(" (%d %d)\n",con_x[i],con_y[i]);
        }
    
        t2=omp_get_wtime();
        if(rank==0)
            printf("time :%lf",t2-t1);
    MPI_Finalize();
    return 0;
}
