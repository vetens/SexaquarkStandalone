#include "../interface/AnalyzerAllSteps.h"


double AnalyzerAllSteps::openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2){
  double opening_angle = TMath::ACos((momentum1.Dot(momentum2))/(pow(momentum1.Mag2()*momentum2.Mag2(),0.5)));
  return opening_angle;
}

double AnalyzerAllSteps::deltaR(double phi1, double eta1, double phi2, double eta2){
	double deltaPhi = reco::deltaPhi(phi1,phi2);
	double deltaEta = eta1-eta2;
	return pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
}


double AnalyzerAllSteps::lxy(TVector3 v1, TVector3 v2){
	double x1 = v1.X();
	double x2 = v2.X();
	double y1 = v1.Y();
	double y2 = v2.Y();
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

double AnalyzerAllSteps::lxyz(TVector3 v1, TVector3 v2){
	double x1 = v1.X();
	double x2 = v2.X();
	double y1 = v1.Y();
	double y2 = v2.Y();
	double z1 = v1.Z();
	double z2 = v2.Z();
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
}



TVector3 AnalyzerAllSteps::PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
   //first move the vector along the line to the starting point of Point_line
   double normalise = sqrt(Vector_along_line.X()*Vector_along_line.X()+Vector_along_line.Y()*Vector_along_line.Y()+Vector_along_line.Z()*Vector_along_line.Z());
   TVector3 n(Vector_along_line.X()/normalise,Vector_along_line.Y()/normalise,Vector_along_line.Z()/normalise);
   TVector3 a = Point_line;
   TVector3 p = Point;

   //see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line (Vector formulation)
   TVector3 vector_PCA = (a-p)-((a-p)*n)*n;
   return vector_PCA ;
}

//return a vector representing the dxy
TVector3 AnalyzerAllSteps::vec_dxy_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){
  //looking at XY, so put the Z component to 0 first
  TVector3 Point_line(Point_line_in.X(),Point_line_in.Y(),0.);
  TVector3 Vector_along_line(Vector_along_line_in.X(), Vector_along_line_in.Y(),0.);
  TVector3 Point(Point_in.X(), Point_in.Y(), 0.);
  
  TVector3 shortest_distance = PCA_line_point(Point_line,  Vector_along_line, Point);
  return shortest_distance;
	
}



double AnalyzerAllSteps::dxy_signed_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){
  TVector3 shortest_distance = vec_dxy_line_point(Point_line_in,Vector_along_line_in,Point_in);
  double dxy_signed_line_point = sqrt(shortest_distance.X()*shortest_distance.X()+shortest_distance.Y()*shortest_distance.Y());

  //looking at XY, so put the Z component to 0 first
  TVector3 Point_line(Point_line_in.X(),Point_line_in.Y(),0.);
  TVector3 Vector_along_line(Vector_along_line_in.X(), Vector_along_line_in.Y(),0.);
  TVector3 Point(Point_in.X(), Point_in.Y(), 0.);

  TVector3 displacement = Point_line - Point; 
  if(displacement*Vector_along_line<0)dxy_signed_line_point = -dxy_signed_line_point;

  return dxy_signed_line_point;
}

double AnalyzerAllSteps::dxyz_signed_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){
  TVector3 shortest_distance = PCA_line_point(Point_line_in,  Vector_along_line_in, Point_in);
  double dxyz_signed_line_point = sqrt(shortest_distance.X()*shortest_distance.X()+shortest_distance.Y()*shortest_distance.Y()+shortest_distance.Z()*shortest_distance.Z());

  TVector3 displacement = Point_line_in - Point_in; 
  if(displacement*Vector_along_line_in<0)dxyz_signed_line_point = -dxyz_signed_line_point;

  return dxyz_signed_line_point;
	
}



double AnalyzerAllSteps::std_dev_lxy(double vx, double vy, double vx_var, double vy_var, double bx_x, double bx_y, double bx_x_var, double bx_y_var){

        double lxy_std_dev_nominator = pow(vx-bx_x,2)*(vx_var+bx_x_var) + pow(vy-bx_y,2)*(vy_var+bx_y_var);
        double lxy_std_dev_denominator = pow(vx-bx_x,2) + pow(vy-bx_y,2);
        double lxy_b_std_dev = sqrt(lxy_std_dev_nominator/lxy_std_dev_denominator);
        return lxy_b_std_dev;

}

//function to return the cos of the angle between the momentum of the particle and it's displacement vector. This is for a V0 particle, so you need the V0 to decay to get it's interaction vertex
double AnalyzerAllSteps::XYpointingAngle(const reco::Candidate  * particle, TVector3 beamspot){
      double angleXY = -2;
      if(particle->numberOfDaughters() == 2){
	      TVector3 decayVertexParticle(particle->daughter(0)->vx(),particle->daughter(0)->vy(),particle->daughter(0)->vz());	 
	      double dx = decayVertexParticle.X()-beamspot.X();
	      double dy = decayVertexParticle.Y()-beamspot.Y();
	      double px = particle->px();
	      double py = particle->py();
	      angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
      }
      return angleXY;
	
}

double AnalyzerAllSteps::CosOpeningsAngle(TVector3 vec1, TVector3 vec2){

  double nom = vec1.X()*vec2.X()+vec1.Y()*vec2.Y()+vec1.Z()*vec2.Z();
  double denom = sqrt(vec1.X()*vec1.X()+vec1.Y()*vec1.Y()+vec1.Z()*vec1.Z())*sqrt(vec2.X()*vec2.X()+vec2.Y()*vec2.Y()+vec2.Z()*vec2.Z());
  return nom/denom;
	
}


double AnalyzerAllSteps::dz_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){
  //looking at Z, so put the XY component to 0 first
//  TVector3 Point_line(0.,0., Point_line_in.Z());
//  TVector3 Vector_along_line(0.,0., Vector_along_line_in.Z());
//  TVector3 Point( 0., 0., Point_in.Z());

//  TVector3 shortest_distance = PCA_line_point(Point_line,  Vector_along_line, Point);
//  return shortest_distance.Z();
  Double_t vz = Point_line_in.Z();
  Double_t vx = Point_line_in.X();
  Double_t vy = Point_line_in.Y();
  Double_t px = Vector_along_line_in.X();
  Double_t py = Vector_along_line_in.Y();
  Double_t pz = Vector_along_line_in.Z(); 
  Double_t pt = sqrt(px*px+py*py);
  //from: https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/TrackBase.h#L764-L767
  return (vz - Point_in.Z()) - ((vx - Point_in.X()) * px + (vy - Point_in.Y()) * py) / pt * pz / pt;
             
}

TVector3 AnalyzerAllSteps::dz_line_point_min(TVector3 Point_line_in, TVector3 Vector_along_line_in, edm::Handle<vector<reco::Vertex>> h_offlinePV){

	TVector3 bestPV;
	double dzmin = 999.;

	for(unsigned int i = 0; i < h_offlinePV->size(); ++i){
		TVector3 PV(h_offlinePV->at(i).x(),h_offlinePV->at(i).y(),h_offlinePV->at(i).z());
		double dz  = AnalyzerAllSteps::dz_line_point(Point_line_in,Vector_along_line_in,PV);
		if(abs(dz) < abs(dzmin)) {dzmin = dz; bestPV =  PV;}
	}

	return bestPV;

}

double AnalyzerAllSteps::sgn(double input){
  double output = 1;
  if(input < 0.) output = -1;
  return output;

}

int AnalyzerAllSteps::getDaughterParticlesTypes(const reco::Candidate * genParticle){
        int pdgIdDaug0 = genParticle->daughter(0)->pdgId();
        int pdgIdDaug1 = genParticle->daughter(1)->pdgId();
        int returnCode = -1;
        if(abs(pdgIdDaug0) == AnalyzerAllSteps::pdgIdPosPion && abs(pdgIdDaug1) == AnalyzerAllSteps::pdgIdPosPion)returnCode = 1;//this is the correct decay mode for Ks to be RECO
        else if((pdgIdDaug0 == AnalyzerAllSteps::pdgIdAntiProton && pdgIdDaug1 == AnalyzerAllSteps::pdgIdPosPion) || (pdgIdDaug1 == AnalyzerAllSteps::pdgIdAntiProton && pdgIdDaug0 == AnalyzerAllSteps::pdgIdPosPion))returnCode = 2;//this is the correct decay mode for an antiLambda to get RECO
        else if((pdgIdDaug0 == AnalyzerAllSteps::pdgIdKs && pdgIdDaug1 == AnalyzerAllSteps::pdgIdAntiLambda) ||(pdgIdDaug1 == AnalyzerAllSteps::pdgIdKs && pdgIdDaug0 == AnalyzerAllSteps::pdgIdAntiLambda)) returnCode = 3;//this is the correct decay mode for an antiS

        return returnCode;

}

int AnalyzerAllSteps::trackQualityAsInt(const reco::Track *track){
    int myquality = -99;
    if(track->quality(reco::TrackBase::undefQuality))myquality = -1;
    if(track->quality(reco::TrackBase::loose))myquality = 0;
    if(track->quality(reco::TrackBase::tight))myquality = 1;
    if(track->quality(reco::TrackBase::highPurity))myquality = 2;
    if(track->quality(reco::TrackBase::confirmed))myquality=3;
    if(track->quality(reco::TrackBase::goodIterative))myquality=4;
    if(track->quality(reco::TrackBase::looseSetWithPV))myquality=5;
    if(track->quality(reco::TrackBase::highPuritySetWithPV))myquality=6;
    if(track->quality(reco::TrackBase::discarded))myquality=7;
    if(track->quality(reco::TrackBase::qualitySize))myquality=8;
    return myquality;
}   
