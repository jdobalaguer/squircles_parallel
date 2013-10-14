function note2(freq, duree, echant)
%%% note(freq, duree) joue une s�rie de notes de fr�quences contenues dans
%%% le tableau freq et de dur�e dans le tableau duree
%%% (en s) (on peut rajouter en 3eme argument la freq d'�chantillonage, par
%%% d�faut 8192Hz

if nargin==2, echant=8192; end
out=[];
for n=1:length(duree);

sery=zeros(1,sum(floor(echant*duree(n))));

ls=0;
sery(ls+1:ls+floor(echant*duree(n)))=sin(2*pi*freq(n)*[1/echant:1/echant:duree(n)]);
ls=ls+floor(echant*duree(n));

out=[out sery];
end

%player=audioplayer(sery, echant);
sound(out,echant);

