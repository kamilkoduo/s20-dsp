// Author (C)  : Samuel GOUGEON - Le Mans, France
// Related Software: Scilab
// License     : CeCILL-B : http://www.cecill.info/licences
// Release     : 2.0
// Date        : 2015-05-24
// Distribution: http://fileexchange.scilab.org/toolboxes/161000
// History     :
//   0.0 - 2010: Creation
//   1.0 - 2011-05-15: first publication on FileExchange
//   2.0 - 2015-05-24: 
//      * usage with structure arrays was undocumented
//      * extension to cell arrays
//
// Reference   : http://bugzilla.scilab.org/show_bug.cgi?id=7293

function R = cshift(M,d)
// CALLING SEQUENCES:
// cshift;          // displays this help
// R = cshift(M, d)
//     circularly shifts by d(i) positions components of M along its #ith dimensions
//
// PARAMETERS:
// M,R : vector, matrix or hypermatrix of any data type. Structure or cell arrays.
// d   : vector of integers. d(i) is the shift to be applied to the M's components 
//        along its #ith dimension.
//
// EXAMPLES:
// // With a matrix or hypermatrix:
// M = resize_matrix([1 2 ; 3 4 ],[3 3 2]);
// M(1:2,1:2,2) = [5 6 ; 7 8 ]
// cshift(M, 1)
// cshift(M, -1)
// cshift(M, [0 1])
// cshift(M, [1 1])
// cshift(M, [0 1 1])
//   
// // With an array of structures:
// s(4,3).b = %t  ;
// s(2,2).r = %pi ;
// s(3,1).t = "test" ;
// cs = cshift(s,[1, -1])
// cs(1,2).b
// cs(3,1).r
// cs(4,3).t
//
// // With a cells array:
// c = cell(4,3);
// c(4,3).entries = %t  ;
// c(2,2).entries = %pi ;
// c(3,1).entries = "test"
// cshift(c,[-1,1])

  Fname = "cshift"
  if argn(2)==0
	head_comments(Fname)
	R = []
	return
  end
  s = size(M)
  R = M
  for i=1:length(d)
    if s(i)>1
      D = pmodulo(d(i),s(i))
      if D~=0
        S = emptystr(1,length(s))+":"
        S(i) = "[s(i)-D+1:s(i) 1:s(i)-D]"
        S = strcat(S,",")
        if typeof(R) ~= "ce"
            execstr("R = R("+S+")")
        else
            execstr("R.entries = R("+S+").entries")
        end
      end
    end
  end
endfunction
