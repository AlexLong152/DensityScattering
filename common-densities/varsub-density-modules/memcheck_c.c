#include <stdio.h>
#include <stdlib.h>
#include "fdefs.h"

#ifdef BGQ
#include <spi/include/kernel/memory.h>
#endif
#ifdef LINUX
#include <memory.h>
#endif
#ifdef MAC
#include <sys/resource.h>
#endif

#ifdef LINUX 
void get_memory_usage_(double *heap_size,double *stack_size,double *heap_avail,double *stack_avail,double *heap_max)
#else
#ifdef BGQ
void get_memory_usage(double *heap_size,double *stack_size,double *heap_avail,double *stack_avail,double *heap_max)
#else
#ifdef MAC 
void get_memory_usage_(double *heap_size,double *stack_size,double *heap_avail,double *stack_avail,double *heap_max)
#else
void get_memory_usage_(double *heap_size,double *stack_size,double *heap_avail,double *stack_avail,double *heap_max)
#endif
#endif
#endif
{
#ifdef LINUX 
  unsigned long memory_size;
  int i,code;
  const char* statm_path = "/proc/self/status";
  const double MB = 1024.0;
  char line[100];
  FILE *f = fopen(statm_path,"r");
  if(!f){
    perror(statm_path);
    abort();
  }
  code=0;
  while(fscanf(f,"%100s",line)!=EOF && code<3)
  {  
   if(strstr(line,"VmPeak")) 
   {
    fscanf(f,"%ld",&memory_size);
    (*heap_max) =memory_size/MB; 
    code++;
   }
   if(strstr(line,"VmSize")) 
   {
    fscanf(f,"%ld",&memory_size);
    (*heap_size) =memory_size/MB; 
    code++;
   }
   if(strstr(line,"VmData")) 
   {
    fscanf(f,"%ld",&memory_size);
    (*stack_size) =memory_size/MB; 
    code++;
   }
  }
  fclose(f);
  
  (*heap_avail) = (double) 0;
  (*stack_avail) =(double) 0; 
#else    
#ifdef BGQ 
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap, heapmax;
  const double MB = 1048576.0;
  
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPMAX, &heapmax);
  
  (*heap_size) =heap/MB;
  (*stack_size) =stack/MB;
  (*heap_avail) =heapavail/MB;
  (*stack_avail) =stackavail/MB;
  (*heap_max) =heapmax/MB;
#else
#ifdef MAC
 const double MB = 1048576.0; 
 struct rusage r_usage;
 
 getrusage(RUSAGE_SELF,&r_usage); 
 (*heap_size) =((double)r_usage.ru_idrss)/MB;
 (*stack_size) =((double)r_usage.ru_isrss)/MB;;
 (*heap_avail) = 0.0;
 (*stack_avail) =0.0;
 (*heap_max) =((double)r_usage.ru_maxrss)/MB; 
#else 
 (*heap_size) =(double) 0;
 (*stack_size) =(double) 0;
 (*heap_avail) =(double) 0;
 (*stack_avail) =(double) 0;
 (*heap_max) =(double) 0;  
#endif 
#endif
#endif
 return;
}
