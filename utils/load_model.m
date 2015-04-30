function model = load_model(organism)
% Load cobra model for the selected organism.
% Options: 'ecoli', 'yeast', 'ecoli_DNApoly', 'yeast_DNApoly'
%
% Author: Daniel Machado, 2013
% u003f 28-apr-15: added DNApoly option

    ECOLI_MODEL = 'models/iAF1260.xml';
    YEAST_MODEL = 'models/iTO977.xml';

    ECOLI_MODEL_DNApoly = 'models/iAF1260_DNApoly.xml';
    YEAST_MODEL_DNApoly = 'models/iTO977_DNApoly.xml';

    switch organism

    	case 'ecoli'
    		model = readCbModel(ECOLI_MODEL);
    	
    	case 'yeast'
    		model = readCbModel(YEAST_MODEL);
    	
    	case 'ecoli_DNApoly'
    		model = readCbModel(ECOLI_MODEL_DNApoly);
    	
    	case 'yeast_DNApoly'
    		model = readCbModel(YEAST_MODEL_DNApoly);
    end   	

    switch organism
        case {'ecoli', 'ecoli_DNApoly'}
            model.lb(model.lb < 0) = -1000;
            model.lb(strcmp('EX_co2(e)', model.rxns)) = 0;
            [model.uptk_rxns, model.secr_rxns] = get_exchange_rxns_ecoli(model);
            model.subtrate_transport_rxn = 'GLCptspp';
            model.subtrate_uptake_rxn = 'EX_glc(e)';

        case {'ecyeastoli', 'yeast_DNApoly'}
            model.c = strcmp('CBIOMASS', model.rxns);
            model.ub(model.ub > 1) = 1000;
            [model.uptk_rxns, model.secr_rxns] = get_exchange_rxns_yeast(model);
            model.subtrate_transport_rxn = 'GAL2_1';
            model.subtrate_uptake_rxn = 'GLCxtI';
            model.lb(strcmp('U214_',model.rxns)) = 0;

            % essential components for anaerobic conditions
            model.ub(strcmp('44DIMZYMSTxtI', model.rxns)) = 1000;
            model.ub(strcmp('C141xtI', model.rxns)) = 1000;
            model.ub(strcmp('C161xtI', model.rxns)) = 1000;
            model.ub(strcmp('C181xtI', model.rxns)) = 1000;
            model.ub(strcmp('ERG572224xtI', model.rxns)) = 1000;
            model.ub(strcmp('LANOSTxtI', model.rxns)) = 1000;
            model.ub(strcmp('ZYMSTxtI', model.rxns)) = 1000;

            % fix problems with lipid composition
            model = addExchangeRxn(model,{'m612'}, 0, 1e-6);
            model = addExchangeRxn(model,{'m709'}, -1000, 0);
    end
end

function [uptk_rxns, secr_rxns] = get_exchange_rxns_ecoli(model)
     uptk_rxns = model.rxns(strncmp('EX', model.rxns, 2) & model.lb < 0 );
     secr_rxns = model.rxns(strncmp('EX', model.rxns, 2) & model.lb == 0 );
end

function [uptk_rxns, secr_rxns] = get_exchange_rxns_yeast(model)
    uptk_rxns = {};
    secr_rxns = {};
    for i = 1:length(model.rxns)
        rxn = model.rxns{i};
        if length(rxn) > 3 && strcmp(rxn(end-2:end), 'xtI')
            uptk_rxns = [uptk_rxns; rxn]; %#ok<AGROW>
        elseif length(rxn) > 3 && strcmp(rxn(end-2:end), 'xtO')
            secr_rxns = [secr_rxns; rxn]; %#ok<AGROW>
        end
    end
end
