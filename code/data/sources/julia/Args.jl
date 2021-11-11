using ArgParse

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--parallel"
    help = "Parallel computation: yes/no"
    arg_type = String
    default = "yes"
    "--nbK"
    help = "Number of sampling points for the 3D Rosenbluth integrals"
    arg_type = Int64
    default = 100
    "--nbAvr"
    help = "Number of sampling points for the orbit-averaging integral"
    arg_type = Int64
    default = 100
    "--q"
    help = "Anisotropy parameter q for the Plummer model."
    arg_type = Float64
    default = 1.0
    ##################################################
    # Options to be removed in GitHub
    ##################################################
    "--gAni"
    help = "Anisotropy parameter gAni for the Isochrone model."
    arg_type = Float64
    default = 0.0
    "--logJrcutoff"
    help = "log10(Jr) lower cut-off"
    arg_type = Float64
    default = -Inf
    "--logLcutoff"
    help = "log10(L) lower cut-off"
    arg_type = Float64
    default = -Inf
end
parsed_args = parse_args(tabargs)
