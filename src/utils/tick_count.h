#ifndef MYTICKCOUNT
#define MYTICKCOUNT

 #include <sys/time.h> /* pour gettimeofday() */
 struct tick_count {
     struct timeval tv_time;
     static tick_count now() {
         tick_count r;
         gettimeofday(& r.tv_time, NULL);
         return r;
     }
     tick_count operator-(tick_count o) const {
         tick_count r;
         /* Perform the carry for the later subtraction by updating y. */
         if (tv_time.tv_usec < o.tv_time.tv_usec) {
             int nsec = (o.tv_time.tv_usec - tv_time.tv_usec) / 1000000 + 1;
             o.tv_time.tv_usec -= 1000000 * nsec;
             o.tv_time.tv_sec += nsec;
         }
         if (tv_time.tv_usec - o.tv_time.tv_usec > 1000000) {
             int nsec = (tv_time.tv_usec - o.tv_time.tv_usec) / 1000000;
             o.tv_time.tv_usec += 1000000 * nsec;
             o.tv_time.tv_sec -= nsec;
         }
         
         /* Compute the time remaining to wait.
          *         tv_usec is certainly posit*ive. */
         r.tv_time.tv_sec = tv_time.tv_sec - o.tv_time.tv_sec;
         r.tv_time.tv_usec = tv_time.tv_usec - o.tv_time.tv_usec;
         
         return r;
     }
     float seconds() { return tv_time.tv_sec + (float)tv_time.tv_usec/1e6; }
 };
 
#endif
 