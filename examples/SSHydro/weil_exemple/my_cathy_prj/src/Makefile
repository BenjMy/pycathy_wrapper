FC90          = gfortran
FC            = gfortran

#FFLAGS        = -O3 -finit-local-zero -Wall
FFLAGS        = -O3
#FFLAGS        = -O2

DEST          = ${HOME}/bin

EXTHDRS       =

HDRS          = CATHY.H IOUNITS.H MB_HGRAPH.H NORMVL.H RIVERNETWORK.H\
		SOILCHAR.H SURFWATER.H TRANSPSURF.H

#LDFLAGS       = -O3 -finit-local-zero -Wall -L/opt/local/lib
#LDFLAGS       = -O3 -L/opt/local/lib
#LDFLAGS       = -O2 -L/opt/local/lib
LDFLAGS       = -O3
#LDFLAGS       = -O2

LIBS          = -lblas -llapack

LINKER        = ${FC}

MAKEFILE      = Makefile

OBJS          = adrstn.o \
		altezze.o \
		anis3d.o \
		area2d.o \
		areas.o \
		assign_dem.o \
		assnew.o \
		asspic.o \
		assut3.o \
		atmbak.o \
		atmnxt.o \
		atmnxt_tra.o \
		atmone.o \
		atmone_tra.o \
		axbsym.o \
		baric_face.o \
		basis3d.o \
		basis6.o \
		bcbak.o \
		bcnew.o \
		bcnxt.o \
		bcnxt_tra.o \
		bcone.o \
		bcone_tra.o \
		bcpic.o \
		bctrn.o \
		bknew.o \
		bkpic.o \
		bkstep.o \
		block_data.o \
		boundir.o \
		calcoli.o \
		cathy_main.o \
		cauchyupd.o \
		cell_nod.o \
		cell_nod_tra.o \
		cflpecnumber.o \
		cfmat.o \
		cfmatn.o \
		cfmatp.o \
		cgsolv.o \
		check.o \
		chknew.o \
		chkpic.o \
		chmaspt.o \
		chmass.o \
		chnew0.o \
		chnew1.o \
		chnew2.o \
		chnew3.o \
		chnew4.o \
		chparm.o \
		chpic0.o \
		chpic1.o \
		chpic2.o \
		chpic3.o \
		chpic4.o \
		chtann.o \
		chtanp.o \
		chvelo.o \
		chvelop.o \
		closio.o \
		conc_bc.o \
		concini_node.o \
		concupd.o \
		concupd_mixing.o \
		conf_fac.o \
		connection_new.o \
		controldiff.o \
		conver.o \
		datin.o \
		datin_tra.o \
		det.o \
		detout.o \
		detoutq.o \
		dirface_update.o \
		dirorneu_face.o \
		dtstat.o \
		durlo.o \
		dxpay.o \
		effnxt.o \
		effone.o \
		eltnod.o \
		eltnod_tra.o \
		etran.o \
		extall.o \
		extcvg.o \
		extr_tvd.o \
		facetype.o \
		fbcdds.o \
		fbcdkr.o \
		fbcdse.o \
		fbckr.o \
		fbcse.o \
		fcdpore.o \
		fcinde.o \
		fhudds.o \
		fhudk2.o \
		fhudk3.o \
		fhudse.o \
		fhukr2.o \
		fhukr3.o \
		fhuse.o \
		flow3d.o \
		fluxmb.o \
		fluxtra.o \
		fvgdds.o \
		fvgdkr.o \
		fvgdse.o \
		fvgkr.o \
		fvgse.o \
		fxvddm.o \
		fxvdkr.o \
		fxvdmc.o \
		fxvkr.o \
		fxvmc.o \
		gen3d.o \
		godunov.o \
		grad.o \
		grdsys.o \
		heapsort.o \
		heron.o \
		hgraph.o \
		icvdwt.o \
		icvhe.o \
		icvhwt.o \
		idmax.o \
		index_dem.o \
		init0.o \
		init0i.o \
		init0r.o \
		init1.o \
		init_sr.o \
		init_sr_tra.o \
		init_surf.o \
		inital.o \
		interp.o \
		interpnew.o \
		interppic.o \
		interppicpt.o \
		interppt.o \
		interpvel.o \
		interpvelpt.o \
		iperplane.o \
		istrato.o \
		kersh.o \
		larson-niklasson.o \
		lhstrn.o \
		locmas.o \
		lsolve.o \
		masbal.o \
		masbil.o \
		mass_adv.o \
		mass_balance.o \
		mass_zone.o \
		masscal.o \
		massmb.o \
		matvett.o \
		mbinit.o \
		mc.o \
		mixing_correction.o \
		mlg_limiter.o \
		moistab.o \
		moistabnew.o \
		moistabpic.o \
		moistabpicpt.o \
		moistabpt.o \
		moistabvel.o \
		moistabvelpt.o \
		neumann.o \
		newton.o \
		newuns.o \
		nod_cell.o \
		nodelt.o \
		nodetotetra.o \
		nodetotetra_tra.o \
		norma2.o \
		normres.o \
		norms.o \
		nudcpt.o \
		nudlocate.o \
		nudnxt.o \
		nudone.o \
		nudpic.o \
		nudwt.o \
		nudwxy.o \
		nudwz.o \
		openio.o \
		ordn.o \
		ordn_v.o \
		pcg.o \
		peatch.o \
		picard.o \
		picuns.o \
		pondupd.o \
		precond.o \
		rast_input_dem.o \
		rast_input_int.o \
		rast_input_lz.o \
		rast_input_real.o \
		rdndbc.o \
		readbc.o \
		recharge.o \
		relax.o \
		relxom.o \
		ressym2.o \
		rhsgrv.o \
		rhsnew.o \
		rhspic.o \
		rhstrn.o \
		riord.o \
		route.o \
		route_tra.o \
		sat_frac.o \
		serbatoio.o \
		sfinit.o \
		sfvbak.o \
		sfvnxt.o \
		sfvone.o \
		sfvrec.o \
		shlnew.o \
		shlpic.o \
		shlsym.o \
		solscal-extended.o \
		sort_elem.o \
		source_trasur.o \
		storcal.o \
		stormb.o \
		strnew.o \
		strpic.o \
		surf_flowtra.o \
		swap_t.o \
		switch.o \
		switch_old.o \
		tetnew.o \
		tetpic.o \
		tim.o \
		time_interpolation.o \
		timnxt.o \
		timupd.o \
		timupdsup.o \
		topia.o \
		tpnodi.o \
		tpnodi2d.o \
		transfer_f3d_surf.o \
		transfer_surf_f3d.o \
                transpmb.o \
		triangoli.o \
		tvd_durlo.o \
		ucdnew.o \
		upd_facetype.o \
		ut3d.o \
		vcopyi.o \
		vcopyr.o \
		vel-rec-utils.o \
		vel3d.o \
		velocity_main.o \
		vettvett.o \
		vnod3d.o \
		volbas.o \
		volfin.o \
		volums.o \
		vtkris3d.o \
		vtkris3dfc.o \
		vtkris3dfcsurf.o \
		weight.o \
		weightneighnode.o \
		wtdepth.o \
		zerorhs.o \
		ztzvec.o

PRINT         = pr

PROGRAM       = cathy_ft

SRCS	      = adrstn.f \
		altezze.f \
		anis3d.f \
		area2d.f \
		areas.f \
		assign_dem.f \
		assnew.f \
		asspic.f \
		assut3.f \
		atmbak.f \
		atmnxt.f \
		atmnxt_tra.f \
		atmone.f \
		atmone_tra.f \
		axbsym.f \
		baric_face.f \
		basis3d.f \
		basis6.f \
		bcbak.f \
		bcnew.f \
		bcnxt.f \
		bcnxt_tra.f \
		bcone.f \
		bcone_tra.f \
		bcpic.f \
		bctrn.f \
		bknew.f \
		bkpic.f \
		bkstep.f \
		block_data.f \
		boundir.f \
		calcoli.f \
		cathy_main.f \
		cauchyupd.f \
		cell_nod.f \
		cell_nod_tra.f \
		cflpecnumber.f \
		cfmat.f \
		cfmatn.f \
		cfmatp.f \
		cgsolv.f \
		check.f \
		chknew.f \
		chkpic.f \
		chmaspt.f \
		chmass.f \
		chnew0.f \
		chnew1.f \
		chnew2.f \
		chnew3.f \
		chnew4.f \
		chparm.f \
		chpic0.f \
		chpic1.f \
		chpic2.f \
		chpic3.f \
		chpic4.f \
		chtann.f \
		chtanp.f \
		chvelo.f \
		chvelop.f \
		closio.f \
		conc_bc.f \
		concini_node.f \
		concupd.f \
		concupd_mixing.f \
		conf_fac.f \
		connection_new.f \
		controldiff.f \
		conver.f \
		datin.f \
		datin_tra.f \
		det.f \
		detout.f \
		detoutq.f \
		dirface_update.f \
		dirorneu_face.f \
		dtstat.f \
		durlo.f \
		dxpay.f \
		effnxt.f \
		effone.f \
		eltnod.f \
		eltnod_tra.f \
		etran.f \
		extall.f \
		extcvg.f \
		extr_tvd.f \
		facetype.f \
		fbcdds.f \
		fbcdkr.f \
		fbcdse.f \
		fbckr.f \
		fbcse.f \
		fcdpore.f \
		fcinde.f \
		fhudds.f \
		fhudk2.f \
		fhudk3.f \
		fhudse.f \
		fhukr2.f \
		fhukr3.f \
		fhuse.f \
		flow3d.f \
		fluxmb.f \
		fluxtra.f \
		fvgdds.f \
		fvgdkr.f \
		fvgdse.f \
		fvgkr.f \
		fvgse.f \
		fxvddm.f \
		fxvdkr.f \
		fxvdmc.f \
		fxvkr.f \
		fxvmc.f \
		gen3d.f \
		godunov.f \
		grad.f \
		grdsys.f \
		heapsort.f \
		heron.f \
		hgraph.f \
		icvdwt.f \
		icvhe.f \
		icvhwt.f \
		idmax.f \
		index_dem.f \
		init0.f \
		init0i.f \
		init0r.f \
		init1.f \
		init_sr.f \
		init_sr_tra.f \
		init_surf.f \
		inital.f \
		interp.f \
		interpnew.f \
		interppic.f \
		interppicpt.f \
		interppt.f \
		interpvel.f \
		interpvelpt.f \
		iperplane.f \
		istrato.f \
		kersh.f \
		larson-niklasson.f \
		lhstrn.f \
		locmas.f \
		lsolve.f \
		masbal.f \
		masbil.f \
		mass_adv.f \
		mass_balance.f \
		mass_zone.f \
		masscal.f \
		massmb.f \
		matvett.f \
		mbinit.f \
		mc.f \
		mixing_correction.f \
		mlg_limiter.f \
		moistab.f \
		moistabnew.f \
		moistabpic.f \
		moistabpicpt.f \
		moistabpt.f \
		moistabvel.f \
		moistabvelpt.f \
		neumann.f \
		newton.f \
		newuns.f \
		nod_cell.f \
		nodelt.f \
		nodetotetra.f \
		nodetotetra_tra.f \
		norma2.f \
		normres.f \
		norms.f \
		nudcpt.f \
		nudlocate.f \
		nudnxt.f \
		nudone.f \
		nudpic.f \
		nudwt.f \
		nudwxy.f \
		nudwz.f \
		openio.f \
		ordn.f \
		ordn_v.f \
		pcg.f \
		peatch.f \
		picard.f \
		picuns.f \
		pondupd.f \
		precond.f \
		rast_input_dem.f \
		rast_input_int.f \
		rast_input_lz.f \
		rast_input_real.f \
		rdndbc.f \
		readbc.f \
		recharge.f \
		relax.f \
		relxom.f \
		ressym2.f \
		rhsgrv.f \
		rhsnew.f \
		rhspic.f \
		rhstrn.f \
		riord.f \
		route.f \
		route_tra.f \
		sat_frac.f \
		serbatoio.f \
		sfinit.f \
		sfvbak.f \
		sfvnxt.f \
		sfvone.f \
		sfvrec.f \
		shlnew.f \
		shlpic.f \
		shlsym.f \
		solscal-extended.f \
		sort_elem.f \
		source_trasur.f \
		storcal.f \
		stormb.f \
		strnew.f \
		strpic.f \
		surf_flowtra.f \
		swap_t.f \
		switch.f \
		switch_old.f \
		tetnew.f \
		tetpic.f \
		tim.f \
		time_interpolation.f \
		timnxt.f \
		timupd.f \
		timupdsup.f \
		topia.f \
		tpnodi.f \
		tpnodi2d.f \
		transfer_f3d_surf.f \
		transfer_surf_f3d.f \
                transpmb.f \
		triangoli.f \
		tvd_durlo.f \
		ucdnew.f \
		upd_facetype.f \
		ut3d.f \
		vcopyi.f \
		vcopyr.f \
		vel-rec-utils.f \
		vel3d.f \
		velocity_main.f \
		vettvett.f \
		vnod3d.f \
		volbas.f \
		volfin.f \
		volums.f \
		vtkris3d.f \
		vtkris3dfc.f \
		vtkris3dfcsurf.f \
		weight.f \
		weightneighnode.f \
		wtdepth.f \
		zerorhs.f \
		ztzvec.f

.SUFFIXES: .prj

.SUFFIXES: .f90

.F.o:
	$(FC90) $(F90FLAGS) -c $<

.f90.o:
	$(FC90) $(F90FLAGS) -c $<

PRJS= $(SRCS:jo=.prj)

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) 
		$(LINKER) -o $(PROGRAM) $(LDFLAGS) $(OBJS) $(LIBS)

check: $(PRJS)
	ftnchek -noverbose -declare $(PRJS) -project -noextern -library > $(PROGRAM).ftn

profile:;       @make "FFLAGS=-pg" "CFLAGS=-pg" "LDFLAGS=-pg" $(PROGRAM)

debug:;         @make "FFLAGS=-g -C -O0 -fbacktrace -fdump-fortran-optimized -fdump-fortran-original" $(PROGRAM)

clean:;		rm -f $(OBJS) $(PROGRAM).cat $(PROGRAM).ftn

depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

print:;		$(PRINT) $(HDRS) $(SRCS)

cat:;		$(CAT) $(HDRS) $(SRCS) > $(PROGRAM).cat

program:        $(PROGRAM)

profile:        $(PROFILE)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@make -f $(MAKEFILE) DEST=$(DEST) install
###
adrstn.o: IOUNITS.H SOILCHAR.H
altezze.o: CATHY.H IOUNITS.H
area2d.o: IOUNITS.H
assign_dem.o: CATHY.H
assnew.o: CATHY.H
asspic.o: CATHY.H
assut3.o: CATHY.H
atmbak.o: CATHY.H 
atmnxt.o: CATHY.H IOUNITS.H  
atmnxt_tra.o: CATHY.H IOUNITS.H
atmone.o: CATHY.H IOUNITS.H SOILCHAR.H 
atmone_tra.o: CATHY.H IOUNITS.H
bcnew.o: CATHY.H
bcpic.o: CATHY.H
bknew.o: CATHY.H
bkpic.o: CATHY.H
bkstep.o: CATHY.H SOILCHAR.H SURFWATER.H IOUNITS.H TRANSPSURF.H
block_data.o: IOUNITS.H
calcoli.o: CATHY.H IOUNITS.H
cathy_main.o: CATHY.H MB_HGRAPH.H SOILCHAR.H SURFWATER.H RIVERNETWORK.H \
	IOUNITS.H TRANSPSURF.H
cauchyupd.o: CATHY.H
cell_nod.o: CATHY.H
check.o: IOUNITS.H
chknew.o: IOUNITS.H
chkpic.o: IOUNITS.H
chmaspt.o: CATHY.H SOILCHAR.H
chmass.o: CATHY.H SOILCHAR.H
chnew0.o: CATHY.H SOILCHAR.H
chnew1.o: CATHY.H SOILCHAR.H
chnew2.o: CATHY.H SOILCHAR.H
chnew3.o: CATHY.H SOILCHAR.H
chnew4.o: CATHY.H SOILCHAR.H
chparm.o: CATHY.H SOILCHAR.H IOUNITS.H
chpic0.o: CATHY.H SOILCHAR.H
chpic1.o: CATHY.H SOILCHAR.H
chpic2.o: CATHY.H SOILCHAR.H
chpic3.o: CATHY.H SOILCHAR.H
chpic4.o: CATHY.H SOILCHAR.H
chtann.o: CATHY.H SOILCHAR.H
chtanp.o: CATHY.H SOILCHAR.H
chvelo.o: CATHY.H SOILCHAR.H
chvelop.o: CATHY.H SOILCHAR.H
closio.o: IOUNITS.H
concini_node.o: IOUNITS.H
concupd.o: IOUNITS.H
connection_new.o: CATHY.H
conver.o: CATHY.H NORMVL.H IOUNITS.H
datin.o: CATHY.H IOUNITS.H SOILCHAR.H SURFWATER.H RIVERNETWORK.H  
datin_tra.o: CATHY.H IOUNITS.H
detout.o: CATHY.H IOUNITS.H
detoutq.o: CATHY.H IOUNITS.H
dirface_update.o: CATHY.H
durlo.o: CATHY.H
effnxt.o: CATHY.H IOUNITS.H
effone.o: CATHY.H IOUNITS.H
eltnod_tra.o: CATHY.H
etran.o: CATHY.H SOILCHAR.H
extall.o: CATHY.H IOUNITS.H
extcvg.o: CATHY.H IOUNITS.H
extr_tvd.o: CATHY.H
facetype.o: CATHY.H
fbcdds.o: CATHY.H SOILCHAR.H
fbcdkr.o: CATHY.H SOILCHAR.H
fbcdse.o: CATHY.H SOILCHAR.H
fbckr.o: CATHY.H SOILCHAR.H
fbcse.o: CATHY.H SOILCHAR.H
fcdpore.o: CATHY.H SOILCHAR.H
fcinde.o: CATHY.H SOILCHAR.H
fhudds.o: CATHY.H SOILCHAR.H
fhudk2.o: CATHY.H SOILCHAR.H
fhudk3.o: CATHY.H SOILCHAR.H
fhudse.o: CATHY.H SOILCHAR.H
fhukr2.o: CATHY.H SOILCHAR.H
fhukr3.o: CATHY.H SOILCHAR.H
fhuse.o: CATHY.H SOILCHAR.H
flow3d.o: CATHY.H NORMVL.H MB_HGRAPH.H IOUNITS.H
fluxmb.o: CATHY.H IOUNITS.H
fluxtra.o: CATHY.H
fvgdds.o: CATHY.H SOILCHAR.H
fvgdkr.o: CATHY.H SOILCHAR.H
fvgdse.o: CATHY.H SOILCHAR.H
fvgkr.o: CATHY.H SOILCHAR.H
fvgse.o: CATHY.H SOILCHAR.H
fxvddm.o: CATHY.H SOILCHAR.H
fxvdkr.o: CATHY.H SOILCHAR.H
fxvdmc.o: CATHY.H SOILCHAR.H
fxvkr.o: CATHY.H SOILCHAR.H
fxvmc.o: CATHY.H SOILCHAR.H
gen3d.o: IOUNITS.H
grdsys.o: IOUNITS.H
hgraph.o: CATHY.H IOUNITS.H SOILCHAR.H
icvdwt.o: IOUNITS.H
icvhe.o: IOUNITS.H
icvhwt.o: IOUNITS.H
index_dem.o: IOUNITS.H
init_sr.o: CATHY.H SURFWATER.H
init_sr_tra.o: CATHY.H TRANSPSURF.H
init_surf.o: CATHY.H SURFWATER.H
inital.o: CATHY.H MB_HGRAPH.H SOILCHAR.H SURFWATER.H RIVERNETWORK.H IOUNITS.H
interp.o: CATHY.H SOILCHAR.H
interpnew.o: CATHY.H SOILCHAR.H
interppic.o: CATHY.H SOILCHAR.H
interppicpt.o: CATHY.H SOILCHAR.H
interppt.o: CATHY.H SOILCHAR.H
interpvel.o: CATHY.H SOILCHAR.H
interpvelpt.o: CATHY.H SOILCHAR.H
larson-niklasson.o: CATHY.H
masbal.o: CATHY.H MB_HGRAPH.H
masbil.o: CATHY.H IOUNITS.H
mbinit.o: IOUNITS.H
mixing_correction.o: CATHY.H
mlg_limiter.o: CATHY.H
neumann.o: CATHY.H
newton.o: CATHY.H IOUNITS.H
nod_cell.o: CATHY.H
nudlocate.o: IOUNITS.H
nudnxt.o: CATHY.H IOUNITS.H
nudone.o: CATHY.H IOUNITS.H
nudpic.o: CATHY.H
openio.o: IOUNITS.H
peatch.o: CATHY.H SOILCHAR.H
picard.o: CATHY.H IOUNITS.H
rast_input_dem.o: CATHY.H
rast_input_int.o: CATHY.H
rast_input_lz.o: CATHY.H
rast_input_real.o: CATHY.H
rhsgrv.o: CATHY.H
route.o: CATHY.H RIVERNETWORK.H SURFWATER.H IOUNITS.H
route_tra.o: CATHY.H RIVERNETWORK.H SURFWATER.H IOUNITS.H TRANSPSURF.H
serbatoio.o: CATHY.H IOUNITS.H
sfinit.o: CATHY.H IOUNITS.H
shlnew.o: CATHY.H
shlpic.o: CATHY.H
shlsym.o: CATHY.H
source_trasur.o: CATHY.H
strnew.o: IOUNITS.H
strpic.o: IOUNITS.H
surf_flowtra.o: CATHY.H SURFWATER.H RIVERNETWORK.H TRANSPSURF.H
switch.o: CATHY.H IOUNITS.H SOILCHAR.H
switch_old.o: CATHY.H SOILCHAR.H
tpnodi.o: CATHY.H IOUNITS.H SOILCHAR.H
tpnodi2d.o: CATHY.H
transfer_f3d_surf.o: CATHY.H
transfer_surf_f3d.o: CATHY.H
transpmb.o: CATHY.H
triangoli.o: CATHY.H
tvd_durlo.o: CATHY.H
upd_facetype.o: CATHY.H
ut3d.o: CATHY.H
vel-rec-utils.o: CATHY.H
vel3d.o: CATHY.H
velocity_main.o: CATHY.H
volbas.o: IOUNITS.H
volfin.o: CATHY.H
wtdepth.o: IOUNITS.H
