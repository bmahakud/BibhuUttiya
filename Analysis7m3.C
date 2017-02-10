#define Analysis7m3_cxx
#include "Analysis7m3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Analysis7m3::Loop()
{//0

TCanvas *c1=new TCanvas("mycanvas1","My Canvas1");
c1->Divide(2,2);

TH1F *h11=new TH1F("Zmass","Z boson mass ",100,40,140);
TH1F *h12=new TH1F("ZPt ","Z Boson Pt",50,0,140);

TH1F *h13=new TH1F("ZPhi","Z boson Phi",50,-5,5);
TH1F *h14=new TH1F("Jet no before deltaR","no of Jet objects before deltaR",40,-2,200);

TCanvas *c2=new TCanvas("mycanvas2","My Canvas2");
c2->Divide(2,2);
TH1F *h21=new TH1F("ZmassS","Z boson mass-91.18 in selected window ",100,-30,30);
TH1F *h22=new TH1F("ZPtselect","ZPt in the selected window ",100,0,140);
TH1F *h23=new TH1F("deltaR","Delta R Jet and electron ",150,-2,10);
TH1F *h24=new TH1F("deltaR after","Delta R Jet and electron after ",150,-2,10);

TCanvas *c3=new TCanvas("mycanvas3","My Canvas3");
c3->Divide(2,2);
TH1F *h31=new TH1F("Jet no after deltaR","no of Jet objects after deltaR",40,-2,200);
TH1F *h32=new TH1F("1st JetPt ","1st Jet Pt for Z boson window events Pt",50,0,140);
TH1F *h33=new TH1F("2nd JetPt ","2nd Jet Pt for Z boson window events Pt",50,0,140);
TH1F *h34=new TH1F("3rd JetPt ","3rd Jet Pt for Z boson window events Pt",50,0,140);

TCanvas *c4=new TCanvas("mycanvas4","My Canvas4");
c4->Divide(2,2);
TH1F *h41=new TH1F("1st JetPt after veto ","1st Jet Pt for Z boson window events after veto",143,-3,140);
TH1F *h42=new TH1F("2nd JetPt after veto ","2nd Jet Pt for Z boson window events after veto",50,0,140);
TH1F *h43=new TH1F("Delta_Phi ","Delta Phi(Z and 1st Jet) after veto",200,-10,10);
TH1F *h44=new TH1F("absDelta_Phi ","fabs(Delta Phi(Z and 1st Jet)) after veto",200,0,10);

TCanvas *c5=new TCanvas("mycanvas5","My Canvas5");
c5->Divide(2,2);
TH1F *h51=new TH1F("1st JetPtaftervetodelataphi ","1st Jet Pt for Z boson window events after veto and delta phi(Z&jet)",50,0,140);
TH1F *h52=new TH1F(" ZPtaftervetodelataphi ","Zboson Pt for Z boson window events after veto and delta phi(Z&jet)",50,0,140);
TH1F *h53=new TH1F("Response","JetPt/zPt",250,-1,6);
TH1F *h54=new TH1F("pfsResponse","pfsJetPt/zPt",250,-1,6);
TCanvas *c6=new TCanvas("mycanvas6","My Canvas6");
c6->Divide(2,2);
TH1F *h61=new TH1F("pfResponse","pfJetPt/zPt",250,-1,6);
TH1F *h62=new TH1F("jptResponse","jptJetPt/zPt",250,-1,6);
TH1F *h63=new TH1F("DeltaRpfcalo","DeltaR between pfjet & calo jet",150,-1,8);

TCanvas *c7=new TCanvas("mycanvas7","My Canvas7");
c7->Divide(2,2);
TH1F *h71=new TH1F("matched calo","matched(with pf) caloJetPt/zPt",250,-1,6);
TH1F *h72=new TH1F("matched pf","matched(with calo) pfJetPt/zPt",250,-1,6);







//   In a ROOT session, you can do:
//      Root > .L Analysis7m3.C
//      Root > Analysis7m3 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

int count=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {//1
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
count=count+1;                     
if(count==100000){                    
break;                             
}

double eE,ePx,ePy,ePz,pE,pPx,pPy,pPz,Pt1,Pt2,zE,zPx,zPy,zPz,zPt,zPhi,zMass;
eE=0;
ePx=0;
ePy=0;
ePz=0;
pE=0;
pPx=0;
pPy=0;
pPz=0;
Pt1=0;
Pt2=0;
zMass=0;
zPhi=0;

int e,p;
e=0;
p=0;

for(int i=0;i<ElePt->size();i++){//zmass,this loop will select GOOD and highest pt electron/positrons from the electrons/positrons
double RelIso= ((EleTrackIso->at(i)+EleCaloIso->at(i))/(ElePt->at(i)));
if(EleCharge->at(i)==-1 &&ElePt->at(i)>Pt1 && ElePt->at(i)>20.0 && fabs(EleEta->at(i)) < 2.5 && EleNChi2->at(i) < 10.0 && RelIso <0.2 && EleIdLoose->at(i)>0){//e-
eE=EleE->at(i);
ePx=ElePx->at(i);
ePy=ElePy->at(i);
ePz=ElePz->at(i);

e=e+1;
Pt1=ElePt->at(i);
}//e-
 if(EleCharge->at(i)==1 && ElePt->at(i)>Pt2 && ElePt->at(i)>20.0 && fabs(EleEta->at(i)) < 2.5 && EleNChi2->at(i) < 10.0 && RelIso <0.2 && EleIdLoose->at(i)>0){//e+
pE=EleE->at(i);
pPx=ElePx->at(i);
pPy=ElePy->at(i);
pPz=ElePz->at(i);

p=p+1;
Pt2=ElePt->at(i);
}//e+
 
}//zmass

if(e>=1 && p>=1){//fill z histo//e>=1,p>=1 means highest pt  e- and e+  are selected in zmass loop

zE=eE+pE;
zPx=ePx+pPx;
zPy=ePy+pPy;
zPz=ePz+pPz;
zPt=sqrt((zPx)*(zPx)+(zPy)*(zPy));
zMass=sqrt((zE*zE)-(zPx*zPx)-(zPy*zPy)-(zPz*zPz));
//define Z phi
zPhi=atan2(zPy,zPx);


h11->Fill(zMass);
h12->Fill(zPt);
h13->Fill(zPhi);






}//fill z histo
int jetn=JetPt->size();
h14->Fill(jetn);

if(fabs(zMass-91.18)<15.0 && jetn>=1){//Z boson events

h21->Fill(zMass-91.18);
h22->Fill(zPt);
//////////////////////////////////////////////////////////////////////calo Jet stuff beggins
std::vector<double> RealCaloJetPt,RealCaloJetEta,RealCaloJetPhi;
RealCaloJetPt.clear();
RealCaloJetPhi.clear();
RealCaloJetEta.clear();

int en=ElePt->size();



double jetEta,jetPhi,eEta,ePhi,DeltaEta,DeltaPhi,DeltaR;
for(int j=0;j<JetPt->size();j++){//jetloop for DeltaR check
jetEta=0;
jetPhi=0;

jetEta=JetEta->at(j);
jetPhi=JetPhi->at(j);
int n=0;
for(int ii=0;ii<ElePt->size();ii++){//elloop
eEta=0;
ePhi=0;
DeltaEta=0;
DeltaPhi=0;
DeltaR=0;
eEta=EleEta->at(ii);
ePhi=ElePhi->at(ii);
DeltaEta=eEta-jetEta;
DeltaPhi=ePhi-jetPhi;
DeltaR=sqrt((DeltaEta)*(DeltaEta)+(DeltaPhi)*(DeltaPhi));
h23->Fill(DeltaR);
if(DeltaR > 0.5){
h24->Fill(DeltaR);
n=n+1;
}

}//elloop


if(n==en){//if all the electrons passed deltaR>0.5 store them as realCalojets



RealCaloJetPt.push_back(JetPt->at(j));
RealCaloJetPhi.push_back(JetPhi->at(j));
RealCaloJetEta.push_back(JetEta->at(j));

}//if all the electrons passed deltaR>0.5 store them as realCalojets



}//jetloop for DeltaR check


//1st jet,2nd jet,3rd jet Ptdistribution
h31->Fill(RealCaloJetPt.size());
//cout<<"real Jet size = "<<RealCaloJetPt.size()<<endl;

for(int ij=0;ij<RealCaloJetPt.size();ij++){//1,2,3jet pt
if(ij==0){//1st jet
h32->Fill(RealCaloJetPt.at(ij));

}//1st jet

if(ij==1){//2nd jet
h33->Fill(RealCaloJetPt.at(ij));

}//2nd jet

if(ij==2){//3rd jet
h34->Fill(RealCaloJetPt.at(ij));

}//3rd jet



}//1,2,3jet pt






if(RealCaloJetPt.size()>=1){//pt balance

if(RealCaloJetPt.size()==1){//1jet event
h41->Fill(RealCaloJetPt.at(0));
h43->Fill(RealCaloJetPhi.at(0)-zPhi);
h44->Fill(fabs(RealCaloJetPhi.at(0)-zPhi));
if(fabs(fabs(RealCaloJetPhi.at(0)-zPhi)-3.14) < 0.3){//delta phi check  jet and Z

h51->Fill(RealCaloJetPt.at(0));
h52->Fill(zPt);
h53->Fill(zPt/(RealCaloJetPt.at(0)));
}//delta phi check jet and Z


}//1jet event
else if(RealCaloJetPt.size()>1){//more than 1 jet
if(RealCaloJetPt.at(1)< 0.2*zPt){//veto
h41->Fill(RealCaloJetPt.at(0));
h42->Fill(RealCaloJetPt.at(1));
h43->Fill(RealCaloJetPhi.at(0)-zPhi);
h44->Fill(fabs(RealCaloJetPhi.at(0)-zPhi));
if(fabs(fabs(RealCaloJetPhi.at(0)-zPhi)-3.14) < 0.3){//delta phi check  jet and Z
h51->Fill(RealCaloJetPt.at(0));
h52->Fill(zPt);
h53->Fill(zPt/(RealCaloJetPt.at(0)));

}//delta phi check jet and Z



}//veto


}//more than 1 jet





}//pt balance

///////////////////////////////////////////////////////////////////////////////////////////////////caloJet stuff ends

////////////////////////////////////////////////////////////////////////////////////////////////////jpt jet stuff beggins
std::vector<double> RealjptJetPt,RealjptJetEta,RealjptJetPhi;
RealjptJetPt.clear();
RealjptJetPhi.clear();
RealjptJetEta.clear();

int jpten=ElePt->size();



double jptjetEta,jptjetPhi,jpteEta,jptePhi,jptDeltaEta,jptDeltaPhi,jptDeltaR;
for(int jjjj=0;jjjj<jptJetPt->size();jjjj++){//pfjetloop for pfDeltaR check
jptjetEta=0;
jptjetPhi=0;

jptjetEta=jptJetEta->at(jjjj);
jptjetPhi=jptJetPhi->at(jjjj);
int jptn=0;
for(int iiii=0;iiii<ElePt->size();iiii++){//pfelloop
jpteEta=0;
jptePhi=0;
jptDeltaEta=0;
jptDeltaPhi=0;
jptDeltaR=0;
jpteEta=EleEta->at(iiii);
jptePhi=ElePhi->at(iiii);
jptDeltaEta=jpteEta-jptjetEta;
jptDeltaPhi=jptePhi-jptjetPhi;
jptDeltaR=sqrt((jptDeltaEta)*(jptDeltaEta)+(jptDeltaPhi)*(jptDeltaPhi));
//h23->Fill(DeltaR);
if(jptDeltaR > 0.5){
//h24->Fill(DeltaR);
jptn=jptn+1;
}

}//elloop


if(jptn==jpten){//if all the electrons passed deltaR>0.5 store them as realCalojets



RealjptJetPt.push_back(jptJetPt->at(jjjj));
RealjptJetPhi.push_back(jptJetPhi->at(jjjj));
RealjptJetEta.push_back(jptJetEta->at(jjjj));

}//if all the electrons passed deltaR>0.5 store them as realCalojets



}//jetloop for DeltaR check


//1st jet,2nd jet,3rd jet Ptdistribution
//h31->Fill(RealCaloJetPt.size());
//cout<<"real Jet size = "<<RealCaloJetPt.size()<<endl;

for(int iijj=0;iijj<RealjptJetPt.size();iijj++){//1,2,3jet pt
if(iijj==0){//1st jet
//h32->Fill(RealCaloJetPt.at(ij));

}//1st jet

if(iijj==1){//2nd jet
//h33->Fill(RealCaloJetPt.at(ij));

}//2nd jet

if(iijj==2){//3rd jet
//h34->Fill(RealCaloJetPt.at(ij));

}//3rd jet



}//1,2,3jet pt






if(RealjptJetPt.size()>=1){//pt balance jpt

if(RealjptJetPt.size()==1){//1jet event
//h41->Fill(RealCaloJetPt.at(0));
//h43->Fill(RealCaloJetPhi.at(0)-zPhi);
//h44->Fill(fabs(RealCaloJetPhi.at(0)-zPhi));
if(fabs(fabs(RealjptJetPhi.at(0)-zPhi)-3.14) < 0.3){//delta phi check  jptjet and Z

//h51->Fill(RealCaloJetPt.at(0));
//h52->Fill(zPt);
//h54->Fill((RealpfJetPt.at(0))/zPt);
h62->Fill(zPt/(RealjptJetPt.at(0)));
}//delta phi check jet and Z


}//1jet event
else if(RealjptJetPt.size()>1){//more than 1 jptjet
if(RealjptJetPt.at(1)< 0.2*zPt){//veto
//h41->Fill(RealpfJetPt.at(0));
//h42->Fill(RealpfJetPt.at(1));
//h43->Fill(RealpfJetPhi.at(0)-zPhi);
//h44->Fill(fabs(RealpfJetPhi.at(0)-zPhi));
if(fabs(fabs(RealjptJetPhi.at(0)-zPhi)-3.14) < 0.3){//delta phi check  jptjet and Z
//h51->Fill(RealpfJetPt.at(0));
//h52->Fill(zPt);
//h54->Fill((RealpfJetPt.at(0))/zPt);
h62->Fill(zPt/(RealjptJetPt.at(0)));
}//delta phi check jptjet and Z



}//veto


}//more than 1 jptjet





}//pt balance jpt






/////////////////////////////////////////////////////////////////////////////////////////////////////jpt jet stuff ends

/////////////////////////////////////////////////////////////////////////////////////////////////////pfjet stuff beggins

std::vector<double> RealpfJetPt,RealpfJetEta,RealpfJetPhi;
RealpfJetPt.clear();
RealpfJetPhi.clear();
RealpfJetEta.clear();

int pfen=ElePt->size();



double pfjetEta,pfjetPhi,pfeEta,pfePhi,pfDeltaEta,pfDeltaPhi,pfDeltaR;
for(int jjj=0;jjj<pfJetPt->size();jjj++){//pfjetloop for pfDeltaR check
pfjetEta=0;
pfjetPhi=0;

pfjetEta=pfJetEta->at(jjj);
pfjetPhi=pfJetPhi->at(jjj);
int pfn=0;
for(int iii=0;iii<ElePt->size();iii++){//pfelloop
pfeEta=0;
pfePhi=0;
pfDeltaEta=0;
pfDeltaPhi=0;
pfDeltaR=0;
pfeEta=EleEta->at(iii);
pfePhi=ElePhi->at(iii);
pfDeltaEta=pfeEta-pfjetEta;
pfDeltaPhi=pfePhi-pfjetPhi;
pfDeltaR=sqrt((pfDeltaEta)*(pfDeltaEta)+(pfDeltaPhi)*(pfDeltaPhi));
//h23->Fill(DeltaR);
if(pfDeltaR > 0.5){
//h24->Fill(DeltaR);
pfn=pfn+1;
}

}//elloop


if(pfn==pfen){//if all the electrons passed deltaR>0.5 store them as realCalojets



RealpfJetPt.push_back(pfJetPt->at(jjj));
RealpfJetPhi.push_back(pfJetPhi->at(jjj));
RealpfJetEta.push_back(pfJetEta->at(jjj));

}//if all the electrons passed deltaR>0.5 store them as realCalojets



}//jetloop for DeltaR check


//1st jet,2nd jet,3rd jet Ptdistribution
//h31->Fill(RealCaloJetPt.size());
//cout<<"real Jet size = "<<RealCaloJetPt.size()<<endl;

for(int iij=0;iij<RealpfJetPt.size();iij++){//1,2,3jet pt
if(iij==0){//1st jet
//h32->Fill(RealCaloJetPt.at(ij));

}//1st jet

if(iij==1){//2nd jet
//h33->Fill(RealCaloJetPt.at(ij));

}//2nd jet

if(iij==2){//3rd jet
//h34->Fill(RealCaloJetPt.at(ij));

}//3rd jet



}//1,2,3jet pt






if(RealpfJetPt.size()>=1){//pt balance pf

if(RealpfJetPt.size()==1){//1jet event
//h41->Fill(RealCaloJetPt.at(0));
//h43->Fill(RealCaloJetPhi.at(0)-zPhi);
//h44->Fill(fabs(RealCaloJetPhi.at(0)-zPhi));
if(fabs(fabs(RealpfJetPhi.at(0)-zPhi)-3.14) < 0.3){//delta phi check  pfjet and Z

//h51->Fill(RealCaloJetPt.at(0));
//h52->Fill(zPt);
//h54->Fill((RealpfJetPt.at(0))/zPt);
h61->Fill(zPt/(RealpfJetPt.at(0)));
}//delta phi check jet and Z


}//1jet event
else if(RealpfJetPt.size()>1){//more than 1 pfjet
if(RealpfJetPt.at(1)< 0.2*zPt){//veto
//h41->Fill(RealpfJetPt.at(0));
//h42->Fill(RealpfJetPt.at(1));
//h43->Fill(RealpfJetPhi.at(0)-zPhi);
//h44->Fill(fabs(RealpfJetPhi.at(0)-zPhi));
if(fabs(fabs(RealpfJetPhi.at(0)-zPhi)-3.14) < 0.3){//delta phi check  pfjet and Z
//h51->Fill(RealpfJetPt.at(0));
//h52->Fill(zPt);
//h54->Fill((RealpfJetPt.at(0))/zPt);
h61->Fill(zPt/(RealpfJetPt.at(0)));
}//delta phi check pfjet and Z



}//veto


}//more than 1 pfjet





}//pt balance pf



///////////////////////////////////////////////////////////////////////////////////////////////////pfJet stuff ends



/////////////////////////////////////////////////
//comparison between calo,jpt and pf starts
/////////////////////////////////////////////


double cje=0;
double cjp=0;
double pje=0;
double pjp=0;
double DR=0;
if(RealCaloJetPt.size()>=1 && RealpfJetPt.size()>=1){//pt balance

if(RealCaloJetPt.size()==1 && RealpfJetPt.size()==1){//1jet event

//if(fabs(fabs(RealCaloJetPhi.at(0)-zPhi)-3.14) < 0.3){//delta phi check  jet and Z
 cje=RealCaloJetEta.at(0);
 cjp=RealCaloJetPhi.at(0);
 pje=RealpfJetEta.at(0);
 pjp=RealpfJetPhi.at(0);
 DR=sqrt(((pjp-cjp)*(pjp-cjp))+((pje-cje)*(pje-cje)));
h63->Fill(sqrt(DR));
if(fabs(fabs(RealCaloJetPhi.at(0)-zPhi)-3.14) < 0.3 && fabs(fabs(RealpfJetPhi.at(0)-zPhi)-3.14) < 0.3 && DR < 0.3){

h71->Fill(zPt/(RealCaloJetPt.at(0)));
h72->Fill(zPt/(RealpfJetPt.at(0)));

}





//}//delta phi check jet and Z


}//1jet event


else if(RealCaloJetPt.size()>1 && RealpfJetPt.size()==1){// >1 calojet & 1 pf
if(RealCaloJetPt.at(1)< 0.2*zPt){//vetocalo
cje=RealCaloJetEta.at(0);
 cjp=RealCaloJetPhi.at(0);
 pje=RealpfJetEta.at(0);
 pjp=RealpfJetPhi.at(0);
 DR=sqrt(((pjp-cjp)*(pjp-cjp))+((pje-cje)*(pje-cje)));
h63->Fill(sqrt(DR));

if(fabs(fabs(RealCaloJetPhi.at(0)-zPhi)-3.14) < 0.3 && fabs(fabs(RealpfJetPhi.at(0)-zPhi)-3.14) < 0.3 && DR < 0.3){

h71->Fill(zPt/(RealCaloJetPt.at(0)));
h72->Fill(zPt/(RealpfJetPt.at(0)));

}


}//vetocalo


}// >1 calojet & 1 pf

else if(RealCaloJetPt.size() ==1 && RealpfJetPt.size()>1){// 1 calojet & >1 pf
if(RealpfJetPt.at(1)< 0.2*zPt){//vetopf

cje=RealCaloJetEta.at(0);
 cjp=RealCaloJetPhi.at(0);
 pje=RealpfJetEta.at(0);
 pjp=RealpfJetPhi.at(0);
 DR=sqrt(((pjp-cjp)*(pjp-cjp))+((pje-cje)*(pje-cje)));
h63->Fill(sqrt(DR));
if(fabs(fabs(RealCaloJetPhi.at(0)-zPhi)-3.14) < 0.3 && fabs(fabs(RealpfJetPhi.at(0)-zPhi)-3.14) < 0.3 && DR < 0.3){

h71->Fill(zPt/(RealCaloJetPt.at(0)));
h72->Fill(zPt/(RealpfJetPt.at(0)));

}

}//vetopf


}// 1 calojet & >1 pf


else if(RealCaloJetPt.size() >1 && RealpfJetPt.size()>1){// >1 calojet & >1 pf
if(RealpfJetPt.at(1)< 0.2*zPt && RealCaloJetPt.at(1)< 0.2*zPt){//vetopfcalo

cje=RealCaloJetEta.at(0);
 cjp=RealCaloJetPhi.at(0);
 pje=RealpfJetEta.at(0);
 pjp=RealpfJetPhi.at(0);
 DR=sqrt(((pjp-cjp)*(pjp-cjp))+((pje-cje)*(pje-cje)));
h63->Fill(sqrt(DR));
if(fabs(fabs(RealCaloJetPhi.at(0)-zPhi)-3.14) < 0.3 && fabs(fabs(RealpfJetPhi.at(0)-zPhi)-3.14) < 0.3 && DR < 0.3){

h71->Fill(zPt/(RealCaloJetPt.at(0)));
h72->Fill(zPt/(RealpfJetPt.at(0)));

}
}//vetopfcalo


}// >1 calojet & >1 pf






}//pt balance














////////////////////////////////////////////////////////
//comparison between calo,jpt and pf ends













}//Z boson events





   }//1

c1->cd(1);
h11->Draw();
c1->cd(2);
h12->Draw();
c1->cd(3);
h13->Draw();
c1->cd(4);
h14->Draw();
c2->cd(1);
h21->Draw();
c2->cd(2);
h22->Draw();
c2->cd(3);
h23->Draw();
c2->cd(4);
h24->Draw();
c3->cd(1);
h31->Draw();
c3->cd(2);
h32->Draw();
c3->cd(3);
h33->Draw();
c3->cd(4);
h34->Draw();
c4->cd(1);
h41->Draw();
c4->cd(2);
h42->Draw();
c4->cd(3);
h43->Draw();
c4->cd(4);
h44->Draw();
c5->cd(1);
h51->Draw();
c5->cd(2);
h52->Draw();
c5->cd(3);
h53->Draw();
//c5->cd(3);
//h54->Draw("SAME");
c6->cd(1);
h61->Draw();
c6->cd(2);
h62->Draw();
c6->cd(3);
h63->Draw();

c7->cd(1);
h71->Draw();
c7->cd(2);
h72->Draw();






}//0



