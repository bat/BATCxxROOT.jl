# This file is a part of BATCxxROOT.jl, licensed under the MIT License (MIT).


function mcmc_model_names(tdir::ATDirectoryInst)
    partreename_rexp = r"^(.*)_parameters$"
    # mcmctreename_rexp = r"^(.*)_mcmc$"
    ks = keys(tdir)
    parname_trees = [String(m[1]) for m in filter(x -> x != nothing, match.(r"^(.*)_parameters$", ks))]
    parname_trees
end


function mcmc_model_names(filename::AbstractString)
    open(TFile, filename, "read") do tfile
        mcmc_model_names(tfile)
    end
end


function mcmc_parameter_names(
    tdir::ATDirectoryInst, model_name::String;
    with_observables::Bool = false,
    with_fixed::Bool = false
)
    tdir_ptr = ROOTFramework.as_pointer(tdir)
    parname_tree_name = "$(model_name)_parameters"

    parnames = icxx"""
        bool with_observables = $with_observables;
        bool with_fixed = $with_fixed;

        TTree* tree = nullptr;
        $tdir_ptr->GetObject($parname_tree_name, tree);

        std::vector<std::string> parnames;

        if (tree->GetBranch("name")) {
            char name[2014];
            bool isparam;
            bool isfixed;
            double fixed_value;
            tree->SetBranchAddress("name", name);
            tree->SetBranchAddress("parameter", &isparam);
            tree->SetBranchAddress("fixed_value", &fixed_value);
            auto nentries = tree->GetEntries();
            for (size_t i = 0; i < nentries; ++i) {
                tree->GetEntry(i);
                // std::cout << "DEBUG: varname = " << name << ", isparam = " << isparam << ", isfixed = " << isfixed << ", fixed_value = " << fixed_value << std::endl;
                if ((with_observables || isparam) && (with_fixed || !isfixed)) {
                    parnames.push_back(name);
                }
            }
        }
        else throw std::runtime_error("Branch name not found");
        parnames;   
    """
end


function read_mcmc_samples(
    filename::String;
    with_observables::Bool = false,
    with_fixed::Bool = false
)
    open(TFile, filename, "read") do tfile
        tfile_ptr = ROOTFramework.as_pointer(tfile)
        model_names = mcmc_model_names(filename)

        model_name = model_names[1]

        parnames_cxx = mcmc_parameter_names(
            tfile, model_name,
            with_observables = with_observables,
            with_fixed = with_fixed
        )
        mcmc_tree_name = "$(model_name)_mcmc"

        params_cxx = icxx"std::vector<double>();"
        logprob_cxx = icxx"std::vector<double>();"
        weight_cxx = icxx"std::vector<int>();"
        chainid_cxx = icxx"std::vector<int32_t>();"
        stepno_cxx = icxx"std::vector<int64_t>();"

        icxx"""
            TFile* tfile = $tfile_ptr;
            auto& parnames = $parnames_cxx;

            auto& params = $params_cxx;
            auto& logprob = $logprob_cxx;
            auto& weight = $weight_cxx;
            auto& chainid = $chainid_cxx;
            auto& stepno = $stepno_cxx;

            TTree* mcmc_tree = nullptr;
            tfile->GetObject($mcmc_tree_name, mcmc_tree);

            size_t nparams = parnames.size();
            std::vector<double> current_params(nparams);
            double current_logprob = 0;
            unsigned current_chain = 0;
            unsigned current_iteration = 0;

            mcmc_tree->SetCacheSize(-1);

            mcmc_tree->SetBranchStatus("*", false);
            for (size_t i = 0; i < nparams; ++i) {
                mcmc_tree->SetBranchStatus(parnames[i].c_str(), true);
                mcmc_tree->SetBranchAddress(parnames[i].c_str(), &current_params[i]);
            }
            mcmc_tree->SetBranchStatus("LogProbability", true);
            mcmc_tree->SetBranchAddress("LogProbability", &current_logprob);
            mcmc_tree->SetBranchStatus("Chain", true);
            mcmc_tree->SetBranchAddress("Chain", &current_chain);
            mcmc_tree->SetBranchStatus("Iteration", true);
            mcmc_tree->SetBranchAddress("Iteration", &current_iteration);

            std::vector<std::vector<double>> last_params;
            std::vector<int> last_weight;
            std::vector<double> last_logprob;
            std::vector<unsigned> last_iteration;

            size_t nsamples = mcmc_tree->GetEntries();
            // std::cout << "DEBUG: nsamples = " << nsamples << std::endl;

            for (size_t i = 0; i < nsamples; ++i) {
                mcmc_tree->GetEntry(i);

                if (last_params.size() < current_chain + 1) {
                    // std::cout << "DEBUG: i = " << i << ", NEW CHAIN: " << current_chain << std::endl;
                    last_params.push_back(std::vector<double>());
                    last_logprob.push_back(0);
                    last_weight.push_back(0);
                    last_iteration.push_back(0);

                    assert(last_params.size() == current_chain + 1);
                    assert(last_logprob.size() == current_chain + 1);
                    assert(last_weight.size() == current_chain + 1);
                    assert(last_iteration.size() == current_chain + 1);

                    last_params[current_chain] = current_params;
                    last_logprob[current_chain] = current_logprob;
                    last_weight[current_chain] = 1;
                    last_iteration[current_chain] = current_iteration;
                } else {
                    // //std::cout << "DEBUG: i = " << i << std::endl;
                    auto& lastpar = last_params[current_chain];
                    // std::cout << "DEBUG: CMP, " << current_params[0] << " " << lastpar[0] << std::endl;
                    if (current_params == lastpar) {
                        assert(last_logprob[current_chain] == current_logprob);
                        last_weight[current_chain] += 1;
                        // std::cout << "DEBUG: INC weight of chain " << current_chain << " to " << last_weight[current_chain] << std::endl;
                    } else {
                        // std::cout << "DEBUG: PUSH params of chain " << current_chain << " with weight " << last_weight[current_chain] << std::endl;

                        params.insert(params.end(), lastpar.begin(), lastpar.end());
                        logprob.push_back(last_logprob[current_chain]);
                        weight.push_back(last_weight[current_chain]);
                        stepno.push_back(last_iteration[current_chain]);
                        chainid.push_back(current_chain + 1);

                        lastpar = current_params;
                        last_logprob[current_chain] = current_logprob;
                        last_weight[current_chain] = 1;
                        last_iteration[current_chain] = current_iteration;
                    }
                }
            }
            for (size_t current_chain = 0; current_chain < last_params.size(); ++current_chain) {
                auto& lastpar = last_params[current_chain];
                if (last_weight[current_chain] > 0) {
                    // std::cout << "DEBUG: final PUSH params of chain " << current_chain << " with weight " << last_weight[current_chain] << std::endl;
                    params.insert(params.end(), lastpar.begin(), lastpar.end());

                    logprob.push_back(last_logprob[current_chain]);
                    weight.push_back(last_weight[current_chain]);
                    stepno.push_back(last_iteration[current_chain]);
                    chainid.push_back(current_chain + 1);
                }
            }
        """

        parnames = map(Symbol, map(String, parnames_cxx))

        nparams = length(parnames)
        nsamples = Int(length(params_cxx) / nparams)

        params = append!(ElasticArray{Float64}(nparams, 0), unsafe_wrap(DenseArray, params_cxx))
        logprob = convert(Vector{Float64}, logprob_cxx)
        weight = convert(Vector{Float64}, weight_cxx)

        @assert size(params, 2) == nsamples
        @assert size(logprob, 1) == nsamples
        @assert size(weight, 1) == nsamples

        smpls = DensitySampleVector(params, logprob, weight)

        chainid = convert(Vector{Int32}, chainid_cxx)
        chaincycle = zeros(chainid)
        stepno = convert(Vector{Int64}, stepno_cxx)
        sampletype = zeros(stepno)

        @assert size(chainid, 1) == nsamples
        @assert size(chaincycle, 1) == nsamples
        @assert size(stepno, 1) == nsamples
        @assert size(sampletype, 1) == nsamples

        smplids = MCMCSampleIDVector(chainid, chaincycle, stepno, sampletype)

        # Free memory:
        icxx"""
            std::vector<double> dummy1;
            std::swap($params_cxx, dummy1);

            std::vector<double> dummy2;
            std::swap($logprob_cxx, dummy2);

            std::vector<int> dummy3;
            std::swap($weight_cxx, dummy3);
        """

        smpls, smplids, parnames
    end
end


#=
struct MCMCSampleID
    chainid::Int32
    chaincycle::Int32
    stepno::Int64
    sampletype::Int64
end
=#
