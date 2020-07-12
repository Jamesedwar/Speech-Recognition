classdef HMM
   
   properties
      Nst
      Nmix
      Qi
      Aij
      Bik
   end
   
   methods (Hidden = true)
      
      % Class constructor
      function H = HMM(Qi,Aij,Bik)
      
      H.Nst  = length(Qi);
      H.Nmix = size(Bik,2);
      H.Qi   = Qi;
      H.Aij  = Aij;
      H.Bik  = Bik;
      
      end
      
      function H0 = getmembers(H)
      
      H0 = struct('Nst',H.Nst,'Nmix',H.Nmix,'Qi',H.Qi,'Aij',H.Aij,'Bik',H.Bik);
      
      end
      
      function [Qi,Aij,Bit,Bitk] = getprob(H,ob)
     
      Qi  = H.Qi;
      Aij = H.Aij;
      
      % normpdf = @(x,m,C) exp(-transpose(x-m)/C*(x-m)/2)/sqrt((2*pi)^length(m)*det(C));
      % normpdf = @(x,m,C) exp(-sum((x-m).^2./diag(C))/2)/sqrt((2*pi)^length(m)*prod(diag(C)));
      % 
      % for i = 1:H.Nst
      %    for t = 1:size(ob,2)
      %       for k = 1:H.Nmix
      %          Bitk(i,t,k) = H.Bik(i,k).w*normpdf(ob(:,t),H.Bik(i,k).m,H.Bik(i,k).C);
      %       end
      %       Bit(i,t) = sum(Bitk(i,t,:));
      %    end
      % end
       
      normpdf = @(x,m,C) exp(-sum(bsxfun(@rdivide,bsxfun(@minus,x,m).^2,diag(C)),1)/2)/sqrt((2*pi)^length(m)*prod(diag(C)));
      for i = 1:H.Nst
         for k = 1:H.Nmix
             Bitk(i,:,k) = H.Bik(i,k).w*normpdf(ob,H.Bik(i,k).m,H.Bik(i,k).C);
         end
      end
      Bit = sum(Bitk,3);
      
      end
      
      function [Qi,Aij,Bit] = getlogprob(H,ob)
      
      [Qi,Aij,Bit] = H.getprob(ob);
      
      Qi  = log(Qi);
      Aij = log(Aij);
      Bit = log(Bit);
      
      end
   end
end
