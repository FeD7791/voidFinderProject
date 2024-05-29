#include "allvars.h"
#include "qsort.h"

void qsort(struct sort *a, int start, int end)
{
  struct sort temp;
  int         left,right,marker;
  double      pivot;

  if (start < end) {

     pivot = a[start].val;
     left = start;
     right = end;

     while (left < right) {
	 while (a[right].val > pivot) {
	     right--;	 
	 }	 

	 while ((left < right) && (a[left].val <= pivot)) {
	     left++;	 
	 }

	 if (left < right) {
	    temp.val = a[left].val;
	    temp.ord = a[left].ord;
            a[left].val = a[right].val;
            a[left].ord = a[right].ord;
            a[right].val = temp.val;	    
            a[right].ord = temp.ord;	    
	 }
     }

     temp.val = a[right].val;
     temp.ord = a[right].ord;
     a[right].val = a[start].val;
     a[right].ord = a[start].ord;
     a[start].val = temp.val;
     a[start].ord = temp.ord;

     marker = right;

     qsort(a, start, marker - 1);
     qsort(a, marker + 1, end);
	  
  }

}	

