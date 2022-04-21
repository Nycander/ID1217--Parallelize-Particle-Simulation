#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>

using namespace std;

#define NUM_THREADS 6

struct thread_data {
   int  thread_id;
   char *message;
};

void *PrintHello(void *threadarg) {
   struct thread_data *my_data;
   my_data = (struct thread_data *) threadarg;
   //sleep(3);

   cout << "Thread ID : " << my_data->thread_id ;
   cout << " Message : " << my_data->message << endl;
   //cout << " Ending thread " ;

   pthread_exit(NULL);
}
void *routine(void *threadarg) {
   printf(" Test ");
}
int main (int argc, char* argv[]) {
   pthread_t threads[NUM_THREADS];
   struct thread_data td[NUM_THREADS];
   int rc;
   int i;
   
   
      cout << "main() start : creating thread[" << 1 << "] next ";
      td[1].thread_id = 1;
      td[1].message = "This is message";
      rc = pthread_create(&threads[1], NULL, PrintHello, (void *)&td[1]);
      
if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
      
      cout << " " ;

      cout << "main() start : creating thread[" << 2 << "] next ";
      td[2].thread_id = 2;
      td[2].message = "This is message";
      rc = pthread_create(&threads[2], NULL, &routine, NULL );
      pthread_join(threads[2], NULL);
      
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }



      cout << "main() start : creating thread[" << 3 << "] next ";
      td[3].thread_id = 3;
      td[3].message = "This is message";
      rc = pthread_create(&threads[3], NULL, PrintHello, (void *)&td[3]);
      
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }



      cout << "main() start : creating thread[" << 4 << "] next ";
      td[4].thread_id = 4;
      td[4].message = "This is message";
      rc = pthread_create(&threads[4], NULL, PrintHello, (void *)&td[4]);
      
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }

      // cout << "creating thread[" << 4 << "] next ";
      // rc = pthread_create(&threads[4], NULL, PrintHello, (void *)4);
      
      // if (rc) {
      //    cout << "Error:unable to create thread," << rc << endl;
      //    exit(-1);
      //    }
   pthread_exit(NULL);
}