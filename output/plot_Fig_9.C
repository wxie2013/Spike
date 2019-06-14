{
    nt->Draw("abs(tw1-tw2)>>h(40, -0.1, 1.1)");
    nt->Draw("abs(w1-w2)>>h2(40, -0.1, 1.1)", "", "same");

    TH1F* r = (TH1F*) h->Clone("r");
    r->Divide(h2);

    c1 = new TCanvas("c1", "", 900, 400);
    c1->Divide(2, 1);

    c1->cd(1);
    c1->SetTitle("distribution");
    h->SetLineColor(2);
    h->Draw();
    h2->Draw("same");
    h->SetTitle("");
    h2->SetTitle("");
    h->GetXaxis()->SetTitle("w1-w2");
    h->GetYaxis()->SetTitle("counts");

    TPaveText *pt = new TPaveText(0.5, 10000, 0.9, 12000);
    pt->AddText("Red: trained");
    pt->AddText("Black: untrained");
    pt->Draw();

    c1->cd(2);
    r->Draw();
    r->SetTitle("");
    r->GetXaxis()->SetTitle("w1-w2");
    r->GetYaxis()->SetTitle("Ratio");

    c1->Update();
}
