{
	gROOT->ProcessLine(".L BetaStudy.C");
	gROOT->ProcessLine("BetaStudy a");
	gROOT->ProcessLine("a.Loop()");
}
