#include<iostream>
#include<math.h>
#include<time.h>
#include<mpi.h>
#include <cstdlib>
using namespace std;
#define MAX 100000 // maximum number of objects to take.

void find_clusters(int n,int k,int *x,int *y,int *kcen){
    struct { 
        double value; 
        int   rank; 
    } kc1[n], o[n],kc[n],d[n];
    int min=999999,max_it=5,it=0,i,j;
    int xcnt=0,cnt=0;
    int cenx[k],ceny[k];
    double euc_distance[k][n];
    int rank,size;
    int root=0,flag=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    srand(time(0));
    for(int i=0;i<n;i++){ //randomly assign ranks to clusters
        kc1[i].rank=rand()%k; //newval //if there is a change is rank
        kc[i].rank=rand()%k; //oldval
        if(i<k)
        {
        cenx[i]=rand()%100;
        ceny[i]=rand()%100;
        }
    }    
    do
    {
        //Calculating the euclidean distances
        for(j=0;j<n;j++) //j is for object number
        {
            d[j].value=(((cenx[rank]-x[j])*(cenx[rank]-x[j]))+((ceny[rank]-y[j])*(ceny[rank]-y[j])));

            d[j].value = sqrt(d[j].value);
            d[j].rank=rank;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //parallel code
        MPI_Reduce( d, kc, n, MPI_DOUBLE_INT, MPI_MINLOC, root,  MPI_COMM_WORLD); //min dist from the centroid to which it belong
        if(rank==root) // update the centroid
        {
        
        int l=0;
        for(j=0;j<n;j++)  //same rank -> same process calculated them.
        {
            if(kc1[j].rank==kc[j].rank)
            {
                l++;
            }
        }
        if(l==n) //l== n means all points on this cluster. so store the rank of the points from kc to kcen
        {
            flag=1;
            for(i=0;i<n;i++)
                kcen[i]=kc[i].rank;
            
        }
        else //if not, put the values of kc to kc1.
        {
            for(j=0;j<n;j++)
                kc1[j].rank=kc[j].rank;
        }

        //Calculating the new centroid co-ordinates
        for(int m=0;m<k;m++) // cluster number
        {
            xcnt=0;
            for(int j=0;j<n;j++) //data objects
            {
                if(kc[j].rank==m)
                {
                    cenx[m]=x[j]+cenx[m];
                    ceny[m]=y[j]+ceny[m];
                    xcnt++;
                }
            }
            //cout<<"\nThe new centroids are: ";
            if(xcnt!=0)
            {
            cenx[m]=(int)(cenx[m]/xcnt);
            ceny[m]=(int)(ceny[m]/xcnt);
            }
            //cout<<cenx[m]<<","<<ceny[m]<<"\n";
        }
    }
        MPI_Bcast(cenx, k, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(ceny, k, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(&flag, 1, MPI_INT, root, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD); // wait for all processes
        it++;

    }while(it<max_it && flag==0);

    //cout<<"it="<<it<<endl;
    if(rank==0)
        for(i=0;i<n;i++) //final rank for each points
            kcen[i]=kc[i].rank;

}
/*int main(int argc,char** argv)
{
//initializing variables
int k,n=100;
int x[n],y[n],kc[n];
int rank, size;
//double d[MAX],euc_distance[MAX][MAX];

//#Clusters
k=3;

//#points
n=100;
srand(time(0));
for(int i=0;i<n;i++)
{
    //x[i]=rand()%100;
    //y[i]=rand()%100;
    x[i]=y[i]=i+1;
    //kc[i]=kc1[i]=0; //initializing as 0
}

MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


find_clusters(n,k,x,y,kc);
//cout<<"hi"<<endl;
if(rank==0)
{
for(int i=0;i<k;i++) //for cluster number
{
    //cout<<"\nCluster "<<i+1<<" centre is :("<<cenx[i]<<","<<ceny[i]<<")";

    cout<<"\nCluster:"<<i+1<<endl;
    int t=0;
    for(int j=0;j<n;j++)
    {
        if(kc[j]==i){ //if the cluster number is equal to i then print the elements i.e. present in the cluster
            cout<<"\n("<<x[j]<<","<<y[j]<<")";
            t+=1;
        }
    }
    cout<<"****"<<t<<endl;
    //MPI_Barrier(MPI_COMM_WORLD);
}

//for(int i=0;i<n;i++)
  //  cout<<x[i]<<" "<<y[i]<<endl;

cout<<"\n";
}
    MPI_Finalize();
return 0;
}*/
