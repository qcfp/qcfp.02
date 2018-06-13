PROJECT=pQCFP


ALLCALCULATORS=  \
calculator-main \
propagatorExciton \
numericalSD 

ALLSUBPROJECTS=  \
asymptoticLF \
averages_orient \
communicator_3rd \
complexv \
constants \
constructor-disorder \
constructor-f-exciton \
constructor-f-exciton-gij \
correlation4points \
dvector3d \
dvectorNd \
dtensor3x3 \
interaction \
eigen \
feinman2sideddiagram1 \
feinman2sideddiagram3 \
interpolationF \
interpolationF2d \
numericalSD \
storage \
toolsFFT \
toolsInterpolate \
toolsInterpolate2d \
toolsIO \
toolsRandom \
toolsCombine \
toolsMatrix \
calculator_redfield \
propagatorODEa \
propagatorM \
propagatorMemory \
propagatorExciton \
propagatorExcitonWF \
calculator_3rd_reduced \
calculator_3rd_wf \
calculator_3rd_secular_cumulant \
calculator_abs_reduced \
calculator_abs_secular_cumulant \
calculator_abs_wf \
calculator-main

all:
	$(foreach val,$(ALLSUBPROJECTS), $(MAKE) -C $(val) &&) $(MAKE) binaries

clean:
	rm -f aaa-bin/* aaa-lib/*



binaries:
	$(foreach val,$(ALLCALCULATORS), cp aaa-bin/z.$(val) aaa-bin//qcfp.$(val).exe &&) pwd

