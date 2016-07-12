#include "IMSRGProfiler.hh"
#include <sys/time.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <omp.h>


map<string, double> IMSRGProfiler::timer;
map<string, int> IMSRGProfiler::counter;
float IMSRGProfiler::start_time = -1;

IMSRGProfiler::IMSRGProfiler()
{
  if (start_time < 0)
  {
    cout << "New profiler..." << endl;
    start_time = omp_get_wtime();
    counter["N_Commutators"] = 0;
    counter["N_Operators"] = 0;
    counter["N_Threads"] = omp_get_max_threads();
    cout << "Initiliation complete! N_threads=" << counter["N_Threads"] << endl;
  }
}
/// Check how much memory is being used.
///
map<string,size_t> IMSRGProfiler::CheckMem()
{
  char cmdstring[100],outbuf[500],buf[100];
  sprintf(cmdstring,"pmap -x %d | tail -1",getpid()); // TODO make this more portable. On OSX, use vmmap. no idea for Windows...
  FILE* output = popen(cmdstring,"r");
  map<string,size_t> s;
  if (output==NULL or fgets(outbuf,500,output) == NULL)
    cout << " <<< IMSRGProfiler::CheckMem():  Problem reading output of pmap (pid = " << getpid() << ")" << endl;
  else
    istringstream(outbuf) >> cmdstring >> buf >> s["Kbytes"] >> s["RSS"] >> s["DIRTY"];
  return s;
}

size_t IMSRGProfiler::MaxMemUsage()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  return (size_t) ru.ru_maxrss;
}

map<string,float> IMSRGProfiler::GetTimes()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  map<string,float> times;
  times["user"] = ru.ru_utime.tv_sec + 1e-6*ru.ru_utime.tv_usec;
  times["system"] = ru.ru_stime.tv_sec + 1e-6*ru.ru_stime.tv_usec;
  times["real"] = omp_get_wtime() - start_time;
  return times;
}

void IMSRGProfiler::PrintTimes()
{
   cout << "Where in the world is this segfault." << endl;
   auto time_tot = GetTimes();
   
   cout << "====================== TIMES (s) ====================" << endl;
   cout.setf(ios::fixed);
   for ( auto it : timer )
   {
     int nfill = (int) (20 * it.second / time_tot["real"]);
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(5) << std::right << it.second;
     cout << " (" << setw(4) << setprecision(1) << 100*it.second / time_tot["real"] << "%) |";
     for (int ifill=0; ifill<nfill; ifill++) cout << "*";
     for (int ifill=nfill; ifill<20; ifill++) cout << " ";
     cout << "|";
     cout  << endl;
     
   }
//   for (auto it : GetTimes())
   for (auto it : time_tot)
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(5) << std::right << it.second  << endl;
}

void IMSRGProfiler::PrintCounters()
{
   cout << "===================== COUNTERS =====================" << endl;
   cout.setf(ios::fixed);
   int N_counters = 0;
   for ( auto jt : counter ) {
      N_counters++;
      cout << "N_Counters=" << N_counters << endl;
      cout << jt.first << endl;
      cout << jt.second << endl;}

   for ( auto it : counter )
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(0) << std::right << it.second  << endl;
   cout << "Got past counters." << endl;
}

void IMSRGProfiler::PrintMemory()
{
   cout << "===================== MEMORY (MB) ==================" << endl;
   for (auto it : CheckMem())
     cout << fixed << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(3) << std::right << it.second/1024. << endl;

   cout << fixed << setw(40) << std::left << "Max Used:  " << setw(12) << setprecision(3) << std::right << MaxMemUsage()/1024.  << endl;
   cout << "Got past Memory." << endl;
}

void IMSRGProfiler::PrintAll()
{
  PrintCounters();
  PrintTimes();
  PrintMemory();
}


