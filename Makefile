ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --libs)

CXXFLAGS      = -O -Wall -fPIC -g -W -Wextra -pedantic
CXXFLAGS     += -isystem $(shell root-config --incdir)
CXXFLAGS     += $(ROOTCFLAGS) 

LIBS	      = $(ROOTLIBS) -lTreePlayer -lTMVA -lRooFitCore -lRooFit -lMinuit -lFoam -lXMLIO

CC            = c++
EXECUTABLE    = run
BINDIR        = ${PWD}/bin/
SRCDIR        = ${PWD}/src/
OBJDIR        = ${PWD}/obj/
#Main program classes
OBJS	      = $(OBJDIR)Main.o $(OBJDIR)Base.o $(OBJDIR)CommonTools.o $(OBJDIR)Fitting.o $(OBJDIR)InternalStorage.o $(OBJDIR)Model.o $(OBJDIR)Pdf_Base.o $(OBJDIR)Pdf_Fit.o $(OBJDIR)Pdf_Gen.o $(OBJDIR)Settings.o $(OBJDIR)Yields.o
#PDF wrappers
OBJS         += $(OBJDIR)Exponential.o $(OBJDIR)myGaussian.o $(OBJDIR)myCrystalBall.o
#Low-mass analytic PDF classes
#DATAOBJS      = $(SRCDIR)RooHILLdini.C #$(SRCDIR)RooHORNSdini.C  $(SRCDIR)RooLITTLEHORNSdini.C 

#-------------------------------------------------------------
#$(EXECUTABLE) :			$(OBJS)	$(OBJDIR)EventDict.o
#				${CC} $(LIBS) -o $(BINDIR)$(EXECUTABLE) $(OBJS) $(OBJDIR)EventDict.o
$(EXECUTABLE) :			$(OBJS)	
				${CC} $(LIBS) -o $(BINDIR)$(EXECUTABLE) $(OBJS)
#-------------------------------------------------------------

#$(OBJDIR)EventDict.o : 		${DATAOBJS}
#				@echo "Generating dictionary ..."
#				@rm -f ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h
#				@rootcint ${SRCDIR}EventDict.C -c ${DATAOBJS}
#				${CC} $(CXXFLAGS) -c ${SRCDIR}EventDict.C -o $(OBJDIR)EventDict.o

$(OBJDIR)Main.o : 		$(SRCDIR)Main.C
				@mkdir -p bin
				${CC} $(CXXFLAGS) -c $(SRCDIR)Main.C -o $(OBJDIR)Main.o

$(OBJDIR)Base.o : 		$(SRCDIR)Base.C $(SRCDIR)Base.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)Base.C -o $(OBJDIR)Base.o

$(OBJDIR)Exponential.o :	$(SRCDIR)Exponential.C $(SRCDIR)Exponential.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)Exponential.C -o $(OBJDIR)Exponential.o

$(OBJDIR)myGaussian.o :	$(SRCDIR)myGaussian.C $(SRCDIR)myGaussian.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)myGaussian.C -o $(OBJDIR)myGaussian.o

$(OBJDIR)myCrystalBall.o :	$(SRCDIR)myCrystalBall.C $(SRCDIR)myCrystalBall.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)myCrystalBall.C -o $(OBJDIR)myCrystalBall.o

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

#-------------------------------------------------------
clean:
		@rm -f $(BINDIR)run $(SRCDIR)/*~ core* $(OBJDIR)/*.o ./*~ ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h
# DO NOT DELETE
