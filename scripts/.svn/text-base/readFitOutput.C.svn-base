void readFitOutput()
{
  TFile f("../output/persistFitOutput.root", "read");

  TMatrixD *CorrelationMatrix = (TMatrixD*)f.Get("CorrelationMatrix");
  TMatrixD *CovarianceMatrix = (TMatrixD*)f.Get("CovarianceMatrix");
  TVectorD *FitParameterValues = (TVectorD*)f.Get("FitParameterValues");
  TVectorD *FitParameterErrors = (TVectorD*)f.Get("FitParameterErrors");
  RooArgList *ListVariablesFloated = (RooArgList*)f.Get("ListVariablesFloated");

  /*
  std::cout << std::endl << "Correlation matrix" << std::endl << std::endl;
  CorrelationMatrix->Print();

  std::cout << std::endl << "Covariance matrix" << std::endl << std::endl;
  CovarianceMatrix->Print();
  */

  std::cout << std::endl << "Fit parameter values" << std::endl << std::endl;
  FitParameterValues->Print();

  std::cout << std::endl << "Fit parameter errors" << std::endl << std::endl;
  FitParameterErrors->Print();

  std::cout << std::endl << "Variables floated" << std::endl << std::endl;
  ListVariablesFloated->Print();

  for (Int_t n_e = 0; n_e < ListVariablesFloated->getSize(); ++n_e)
    {
      std::cout << ListVariablesFloated->at(n_e)->GetName() << " : "
                << (*FitParameterValues)(n_e) << " +/- " << (*FitParameterErrors)(n_e)
                << std::endl;
    }

}
