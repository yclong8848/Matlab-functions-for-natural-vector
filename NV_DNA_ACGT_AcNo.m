function [NV] = NV_DNA_ACGT_AcNo(AcNo) %The computer must be online, example: AcNo='XM_012746469.1'

point = strfind(AcNo, '.');

if  point > 1
   
    LOCUS=AcNo(1:(point-1));

else  LOCUS=AcNo;
    
end

sequ = getgenbank(LOCUS, 'SequenceOnly', true);

S = nt2int(sequ);
% Transform the sequence into integers (Matlab):
% Nucleotide	                         Code	Integer
% Adenosine 	                          A	       1
% Cytidine   	                          C	       2
% Guanine 	                              G	       3
% Thymidine 	                          T	       4
% Uridine (if 'Alphabet' set to 'RNA') 	  U	       4
% Purine (A or G) 	                      R	       5
% Pyrimidine (T or C) 	                  Y	       6
% Keto (G or T) 	                      K	       7
% Amino (A or C) 	                      M	       8
% Strong interaction (3 H bonds) (G or C) S	       9
% Weak interaction (2 H bonds) (A or T)   W	      10
% Not A (C or G or T)	                  B	      11
% Not C (A or G or T)	                  D	      12
% Not G (A or C or T)	                  H	      13
% Not T or U (A or C or G)	              V	      14
% Any nucleotide (A or C or G or T or U)  N	      15
% Gap of indeterminate length	          -	      16 (ignore here)
% Unknown (any character not in table)	  *	       0 (ignore here)


na=0;
Ta=0;
ta=[];  % Weighted values for all A in the sequence

nc=0;
Tc=0;
tc=[];  % Weighted values for all C in the sequence

ng=0;
Tg=0;
tg=[];  % Weighted values for all G in the sequence

nt=0;
Tt=0;
tt=[];  % Weighted values for all T in the sequence
% Initial values

nn= length(S); 


for i= 1: nn
    
    if S(i)==1;
        na=na+1;
        Ta=Ta+i;    %Treat the null base before the 1st base as the origin %
        ta(i)=1;    
        tc(i)=0;
        tg(i)=0;
        tt(i)=0;
    elseif S(i)==2;
        nc=nc+1;
        Tc=Tc+i;   
        ta(i)=0;    
        tc(i)=1;
        tg(i)=0;
        tt(i)=0;
    elseif S(i)==3;
        ng=ng+1;
        Tg=Tg+i;       
        ta(i)=0;    
        tc(i)=0;
        tg(i)=1;
        tt(i)=0;
    elseif S(i)==4;
        nt=nt+1;
        Tt=Tt+i;   
        ta(i)=0;    
        tc(i)=0;
        tg(i)=0;
        tt(i)=1;
        
        
    elseif S(i)==5; 
        na=na+(1/2);
        Ta=Ta+ i*(1/2);      
        ng=ng+(1/2);
        Tg=Tg+ i*(1/2);   
        ta(i)=(1/2);    
        tc(i)=0;
        tg(i)=(1/2);
        tt(i)=0;
    elseif S(i)==6; 
        nt=nt+(1/2);
        Tt=Tt+ i*(1/2);      
        nc=nc+(1/2);
        Tc=Tc+ i*(1/2);    
        ta(i)=0;    
        tc(i)=(1/2);
        tg(i)=0;
        tt(i)=(1/2);   
    elseif S(i)==7; 
        ng=ng+(1/2);
        Tg=Tg+ i*(1/2);   
        nt=nt+(1/2);
        Tt=Tt+ i*(1/2);   
        ta(i)=0;    
        tc(i)=0;
        tg(i)=(1/2);
        tt(i)=(1/2);   
    elseif S(i)==8; 
        na=na+(1/2);
        Ta=Ta+ i*(1/2);      
        nc=nc+(1/2);
        Tc=Tc+ i*(1/2);    
        ta(i)=(1/2);    
        tc(i)=(1/2);
        tg(i)=0;
        tt(i)=0;   
    elseif S(i)==9; 
        ng=ng+(1/2);
        Tg=Tg+ i*(1/2);     
        nc=nc+(1/2);
        Tc=Tc+ i*(1/2);   
        ta(i)=0;    
        tc(i)=(1/2);
        tg(i)=(1/2);
        tt(i)=0;   
    elseif S(i)==10; 
        na=na+(1/2);
        Ta=Ta+ i*(1/2);    
        nt=nt+(1/2);
        Tt=Tt+ i*(1/2);     
        ta(i)=(1/2);    
        tc(i)=0;
        tg(i)=0;
        tt(i)=(1/2);    
        
        
    elseif S(i)==11; 
        nc=nc+(1/3);
        Tc=Tc+ i*(1/3);   
        ng=ng+(1/3);
        Tg=Tg+ i*(1/3);   
        nt=nt+(1/3);
        Tt=Tt+ i*(1/3);   
        ta(i)=0;    
        tc(i)=(1/2);
        tg(i)=(1/2);
        tt(i)=(1/2);    
    elseif S(i)==12; 
        na=na+(1/3);
        Ta=Ta+ i*(1/3);     
        ng=ng+(1/3);
        Tg=Tg+ i*(1/3);   
        nt=nt+(1/3);
        Tt=Tt+ i*(1/3);   
        ta(i)=(1/3);    
        tc(i)=0;
        tg(i)=(1/3);
        tt(i)=(1/3);         
    elseif S(i)==13; 
        na=na+(1/3);
        Ta=Ta+ i*(1/3);   
        nc=nc+(1/3);
        Tc=Tc+ i*(1/3);   
        nt=nt+(1/3);
        Tt=Tt+ i*(1/3);    
        ta(i)=(1/3);    
        tc(i)=(1/3);
        tg(i)=0;
        tt(i)=(1/3);         
    elseif S(i)==14; 
        na=na+(1/3);
        Ta=Ta+ i*(1/3);     
        nc=nc+(1/3);
        Tc=Tc+ i*(1/3);     
        ng=ng+(1/3);
        Tg=Tg+ i*(1/3);      
        ta(i)=(1/3);    
        tc(i)=(1/3);
        tg(i)=(1/3);
        tt(i)=0;      
        
        
    elseif S(i)==15; 
        na=na+(1/4);
        Ta=Ta+ i*(1/4);     
        nc=nc+(1/4);
        Tc=Tc+ i*(1/4);    
        ng=ng+(1/4);
        Tg=Tg+ i*(1/4);        
        nt=nt+(1/4);
        Tt=Tt+ i*(1/4);   
        ta(i)=(1/4);    
        tc(i)=(1/4);
        tg(i)=(1/4);
        tt(i)=(1/4);     
        
    end
    
end

% Here are the mu values of A, C, G, T.

ua=Ta/na;                    
uc=Tc/nc;
ug=Tg/ng;
ut=Tt/nt;

% Here are the variances of A, C, G, T for natural vector.

j=2;                         % Choose j order variance for natural vector.
    
Da=0;
Dc=0;
Dg=0;
Dt=0;
% Initial values

for i=1:nn 

    Da=Da+ ((i-ua)^j)*ta(i)/((na^(j-1))*(nn^(j-1)));
    
    Dc=Dc+ ((i-uc)^j)*tc(i)/((nc^(j-1))*(nn^(j-1)));
    
    Dg=Dg+ ((i-ug)^j)*tg(i)/((ng^(j-1))*(nn^(j-1)));
    
    Dt=Dt+ ((i-ut)^j)*tt(i)/((nt^(j-1))*(nn^(j-1)));
    

    NV =[na,nc,ng,nt, ua,uc,ug,ut, Da,Dc,Dg,Dt];
    
     
end
    
    
    
    