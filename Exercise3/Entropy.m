function H = Entropy( TypeEntropy, p, q )
%
%	H = Entropy( TypeEntropy, p, q  )   {or H = Entropy( p ) }
%
%	Generates the Entropy calculation for a given vector or matrix.  
%
%	If TypeEntropy is not given, normal Entropy calulation based apon
%		H(X) = - sum( p(x) * log( p(x) ) )
%
%	For TypeEntropy = 'J'
%		Joint Entropy 
%			-> H(X,Y) = - double sum( p(x,y) * log( p(x,y) ) )
%
%	For TypeEntropy = 'C'
%		Conditional Entropy 
%			-> H(X|Y) = H(X,Y) - H(Y)
%					Input is a matrix with X along the top and Y along the side
%					(To get H(Y|X), simply send the Transpose of X to function)
%
%	For TypeEntropy = 'M'
%		Marginal Entropy 
%			-> H(X) = - sum( p(x) * log( p(x) ) )
%					Input is a matrix with X along the top and Y along the side
%					(To get H(Y), simply send the Transpose of p(X) to function)
%
%	For TypeEntropy = 'R'
%		Relative Entropy (Kullback Leibler Distance)
%			-> H(X,Y) = D = sum( p(p>0) .* log2( p(p>0)./q(q>0) ));
%
%	For TypeEntropy = 'I'
%		Mutual Information
%			-> I(X;Y) = - double sum( p(x,y) * log( p(x,y)/(p(x)p(y)) ) )
%			-> I(X;Y) =  H(X) - H(X|Y)
%					Input is a matrix with X along the top and Y along the side
%					(To get I(Y;X), simply send the Transpose of p(X) to function)
%
%	========================================================
%	Copyright 1999, Philip M. Hanna, Wright State University
%	For use in Information Theory - EE 740, with permission
%	of Fred Garber, Ph.D.
%	========================================================
%
if nargin == 0,
   help Entropy
   H = [];
   
elseif nargin == 1,
   p = TypeEntropy;
   H = -sum( p(p>0) .* log2( p(p>0) ))';
   
else
	switch lower(TypeEntropy)
   case 'j', 
      H = -sum( sum (p(p>0) .* log2( p(p>0) )));
      
   case 'c',
%
%	Conditional Entropy
%			->  H(X|Y) = H(X,Y) - H(Y)
		H = Entropy('j',p) - Entropy('m',p');      
      
   case 'm'
%	Marginal Entropy
      MarginalX = sum(p);
		H = -sum( MarginalX(MarginalX>0) .* log2( MarginalX(MarginalX>0) ));
      
   case 'r'
%	Kullback Leibler Distance      
%  D = sum( p(p>0) .* log2( p(p>0)./q(q>0) )); 
      p = p' + eps;
   	q = q' + eps;
   	H = sum( p .* log2( p./q ));
      
   case 'i'
%		Mutual Information
%			-> I(X;Y) = - double sum( p(x,y) * log( p(x,y)/(p(x)p(y)) ) )
%			-> I(X;Y) =  H(X) - H(X|Y)
		H   = Entropy('m',p) - Entropy('c',p);
      
   otherwise, 
      disp('Unknown method.')
   end
end