void new2(TString fname1, TString fname2)
{
	TFile f1(fname1);
	TTree* t1 = (TTree*) f1.Get("EDepSimEvents");
	TFile f2(fname1);
	TTree* t2 = (TTree*) f2.Get("EDepSimEvents");
	
	TG4Event* ev1 = new TG4Event;
	TG4Event* ev2 = new TG4Event;
	
    t1->SetBranchAddress("Event",&ev1);
    t2->SetBranchAddress("Event",&ev2);
	
	for(int i = 0; i < t1->GetEntries(); i++)
	{
		t1->GetEntry(i);
		t2->GetEntry(i);
		
		for(int j = 0; j < ev1->SegmentDetectors["EMCalSci"].size(); j++)
		{
				TG4HitSegment hit1 = ev1->SegmentDetectors["EMCalSci"].at(j);
				TG4HitSegment hit2 = ev2->SegmentDetectors["EMCalSci"].at(j);
				
				if(hit1.GetStart().X() != hit2.GetStart().X() ||
				hit1.GetStart().Y() != hit2.GetStart().Y() ||
				hit1.GetStart().Z() != hit2.GetStart().Z())
				{
					std::cout << hit1.GetStart().X() << " " << hit1.GetStart().Y() << " " << hit1.GetStart().Z() << std::endl;
				}
		}	
	}
}