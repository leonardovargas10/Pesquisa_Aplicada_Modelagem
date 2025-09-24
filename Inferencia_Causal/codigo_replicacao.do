


* Definir a pasta de trabalho
cd "SUA PASTA AQUI"


* Replicacao of Dranove, Ody, Starc (2021)
clear all
use if sample1==1 using "dataset", clear


* Criando grupos de tratamento
gen gvar_aux = yq if L0 == 1
bysort state_r: egen gvar_CS = mean(gvar_aux)
replace gvar_CS = 0 if gvar_CS == .


* lista de variaveis de interesse
global dvlist lnspend_pc lnpriceperpresc lnpresc_pc lnquantity_pcp genpen_presc ///
       geneff_presc genacc_presc genacc_sim s_MCO s_offset


* loop para estimar o DiD para todas as variaveis
foreach depvar of global dvlist { 	

	* CSDiD never treated
	csdid2 `depvar' exp [w=wgt], gvar(gvar_CS) tvar(yq) ivar(state_r) agg(event) long2
	estat event, window(-9 9) estore(cs_nt_`depvar')
	
	* CSDiD not yet treated
	csdid2 `depvar' exp [w=wgt], gvar(gvar_CS) tvar(yq) ivar(state_r) agg(event) long2 notyet
	estat event, window(-9 9) estore(cs_nyt_`depvar')
}




* Graficos do Event-Study
foreach depvar in $dvlist {
    event_plot cs_nt_`depvar' cs_nyt_`depvar', ///
        stub_lag(tp# tp#) stub_lead(tm# tm#) ///
		together perturb(-0.3(0.1)0.3) noautolegend ///
		plottype(scatter) ciplottype(rcap) ///
            lag_opt1(msymbol(lgx) msize(1.3) mlwidth(0.5) ///
			color(green)) lag_ci_opt1(color(green) lw(0.5)) ///
            lag_opt2(msymbol(Dh) msize(1.3) mlwidth(0.5) ///
			color(orange)) lag_ci_opt2(color(orange) lw(0.5)) ///
        graph_opt(xtitle("") ytitle("") name(`depvar', replace) ///
			subtitle(" `:variable label `depvar''") xlabel(-9(1)9) ///
			legend(order(1 "CSDiD Never Treated" 3 "CSDiD Not Yet Treated") ///
			size(medium) pos(6) rows(1) region(style(none))) ///
			xline(-0.5, lc(gs8) lp(dash)) yline(0, lc(gs8) lp(dash)))

graph export "replication_results/id_51 - `:variable label `depvar'' stata.pdf", replace
}
