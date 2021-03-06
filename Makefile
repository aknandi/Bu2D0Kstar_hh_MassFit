ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --libs)

CXXFLAGS      = -O -Wall -fPIC -g -W -Wextra -pedantic
CXXFLAGS     += -isystem $(shell root-config --incdir)
CXXFLAGS     += $(ROOTCFLAGS) 

LIBS	      = $(ROOTLIBS) -lTreePlayer -lTMVA -lRooFitCore -lRooFit -lRooStats -lMinuit -lFoam -lXMLIO

CC            = c++
EXECUTABLE    = run
BINDIR        = ${PWD}/bin/
SRCDIR        = ${PWD}/src/
OBJDIR        = ${PWD}/obj/
#Main program classes
OBJS	      = $(OBJDIR)Main.o $(OBJDIR)Base.o $(OBJDIR)CommonTools.o $(OBJDIR)Fitting.o $(OBJDIR)InternalStorage.o $(OBJDIR)Model.o $(OBJDIR)Pdf_Base.o $(OBJDIR)Pdf_Fit.o $(OBJDIR)Pdf_Gen.o $(OBJDIR)Settings.o $(OBJDIR)Yields.o
#PDF wrappers
OBJS         += $(OBJDIR)Exponential.o $(OBJDIR)myGaussian.o $(OBJDIR)DoubleGaussian.o $(OBJDIR)DoubleCrystalBall.o $(OBJDIR)DoubleJohnson.o $(OBJDIR)PartRecoDstKst.o $(OBJDIR)PartRecoShapes.o $(OBJDIR)myCruijff.o
#Low-mass analytic PDF classes
#DATAOBJSCINT  = RooHORNSdini.C RooHILLdini.C RooLITTLEHORNSdini.C 
DATAOBJS      = $(SRCDIR)RooHORNSdini.C $(SRCDIR)RooHILLdini.C $(SRCDIR)RooLITTLEHORNSdini.C $(SRCDIR)RooJohnsonSU.C $(SRCDIR)RooCruijff.C 
#-------------------------------------------------------------
$(EXECUTABLE) :			$(OBJS)	$(OBJDIR)EventDict.o
				${CC} $(LIBS) -o $(BINDIR)$(EXECUTABLE) $(OBJS) $(OBJDIR)EventDict.o
#$(EXECUTABLE) :			$(OBJS)	
#				${CC} $(LIBS) -o $(BINDIR)$(EXECUTABLE) $(OBJS)
#-------------------------------------------------------------

$(OBJDIR)EventDict.o : 		${DATAOBJS}
				@echo "Generating dictionary ..."
				@rm -f ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h
				@rootcint ${SRCDIR}EventDict.C -c ${DATAOBJS}
				${CC} $(CXXFLAGS) -I ${PWD} -c ${SRCDIR}EventDict.C -o $(OBJDIR)EventDict.o

$(OBJDIR)Main.o : 		$(SRCDIR)Main.C
				@mkdir -p bin
				${CC} $(CXXFLAGS) -c $(SRCDIR)Main.C -o $(OBJDIR)Main.o

$(OBJDIR)Base.o : 		$(SRCDIR)Base.C $(SRCDIR)Base.h
						${CC} $(CXXFLAGS) -c $(SRCDIR)Base.C -o $(OBJDIR)Base.o

$(OBJDIR)Exponential.o :	$(SRCDIR)Exponential.C $(SRCDIR)Exponential.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)Exponential.C -o $(OBJDIR)Exponential.o

$(OBJDIR)myGaussian.o :	$(SRCDIR)myGaussian.C $(SRCDIR)myGaussian.h
						${CC} $(CXXFLAGS) -c $(SRCDIR)myGaussian.C -o $(OBJDIR)myGaussian.o

$(OBJDIR)DoubleGaussian.o :	$(SRCDIR)DoubleGaussian.C $(SRCDIR)DoubleGaussian.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)DoubleGaussian.C -o $(OBJDIR)DoubleGaussian.o
							
$(OBJDIR)DoubleCrystalBall.o :	$(SRCDIR)DoubleCrystalBall.C $(SRCDIR)DoubleCrystalBall.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)DoubleCrystalBall.C -o $(OBJDIR)DoubleCrystalBall.o

$(OBJDIR)DoubleJohnson.o :	$(SRCDIR)DoubleJohnson.C $(SRCDIR)DoubleJohnson.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)DoubleJohnson.C -o $(OBJDIR)DoubleJohnson.o
							
$(OBJDIR)CommonTools.o : 	$(SRCDIR)CommonTools.C $(SRCDIR)CommonTools.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)CommonTools.C -o $(OBJDIR)CommonTools.o

$(OBJDIR)Fitting.o : 		$(SRCDIR)Fitting.C $(SRCDIR)Fitting.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)Fitting.C -o $(OBJDIR)Fitting.o

$(OBJDIR)InternalStorage.o : 	$(SRCDIR)InternalStorage.C $(SRCDIR)InternalStorage.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)InternalStorage.C -o $(OBJDIR)InternalStorage.o

$(OBJDIR)Model.o : 		$(SRCDIR)Model.C $(SRCDIR)Model.h
						${CC} $(CXXFLAGS) -c $(SRCDIR)Model.C -o $(OBJDIR)Model.o

$(OBJDIR)Pdf_Base.o : 		$(SRCDIR)Pdf_Base.C $(SRCDIR)Pdf_Base.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)Pdf_Base.C -o $(OBJDIR)Pdf_Base.o

$(OBJDIR)Pdf_Fit.o : 		$(SRCDIR)Pdf_Fit.C $(SRCDIR)Pdf_Fit.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)Pdf_Fit.C -o $(OBJDIR)Pdf_Fit.o

$(OBJDIR)Pdf_Gen.o : 		$(SRCDIR)Pdf_Gen.C $(SRCDIR)Pdf_Gen.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)Pdf_Gen.C -o $(OBJDIR)Pdf_Gen.o

$(OBJDIR)Yields.o : 		$(SRCDIR)Yields.C $(SRCDIR)Yields.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)Yields.C -o $(OBJDIR)Yields.o

$(OBJDIR)Settings.o : 		$(SRCDIR)Settings.C $(SRCDIR)Settings.h
							${CC} $(CXXFLAGS) -c $(SRCDIR)Settings.C -o $(OBJDIR)Settings.o
				
$(OBJDIR)PartRecoShapes.o :     $(SRCDIR)PartRecoShapes.C $(SRCDIR)PartRecoShapes.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)PartRecoShapes.C -o $(OBJDIR)PartRecoShapes.o
                                
$(OBJDIR)PartRecoDstKst.o :             $(SRCDIR)PartRecoDstKst.C $(SRCDIR)PartRecoDstKst.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)PartRecoDstKst.C -o $(OBJDIR)PartRecoDstKst.o

$(OBJDIR)RooHORNSdini.o :       $(SRCDIR)RooHORNSdini.C $(SRCDIR)RooHORNSdini.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)RooHORNSdini.C -o $(OBJDIR)RooHORNSdini.o
								
$(OBJDIR)RooHILLdini.o :        $(SRCDIR)RooHILLdini.C $(SRCDIR)RooHILLdini.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)RooHILLdini.C -o $(OBJDIR)RooHILLdini.o

$(OBJDIR)RooLITTLEHORNSdini.o :       $(SRCDIR)RooLITTLEHORNSdini.C $(SRCDIR)RooLITTLEHORNSdini.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)RooLITTLEHORNSdini.C -o $(OBJDIR)RooLITTLEHORNSdini.o
								
$(OBJDIR)RooJohnsonSU.o :       $(SRCDIR)RooJohnsonSU.C $(SRCDIR)RooJohnsonSU.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)RooJohnsonSU.C -o $(OBJDIR)RooJohnsonSU.o	
								
$(OBJDIR)RooCruijff.o :     $(SRCDIR)RooCruijff.C $(SRCDIR)RooCruijff.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)RooCruijff.C -o $(OBJDIR)RooCruijff.o						

$(OBJDIR)myCruijff.o :     $(SRCDIR)myCruijff.C $(SRCDIR)myCruijff.h
								${CC} $(CXXFLAGS) -c $(SRCDIR)myCruijff.C -o $(OBJDIR)myCruijff.o	
#-------------------------------------------------------
clean:
		@rm -f $(BINDIR)run $(SRCDIR)/*~ core* $(OBJDIR)/*.o ./*~ ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h
# DO NOT DELETE
