from optparse import OptionParser

def initParserCommandLine():
    parser = OptionParser()
    parser.add_option("-c", "--cores", dest="cores", default = 10,
                  help="Count of core")

    parser.add_option("-s", "--sites", default = 611, dest="sites",
                  help="Count of sites")

    parser.add_option("-d", "--diffuse", default = 0.0001, dest="diffuse",
                  help="Parametr of diffuse")

    parser.add_option("-q", "--trans_start", default = 0.002069, dest="trans_start",
                  help="Parametr of transcription start")

    parser.add_option("-w", "--kni_degrad", default = 0.000009627, dest="kni_degrad",
                  help="Parametr of kni degradation")

    parser.add_option("-e", "--mrna_degrad", default = 0.000019254, dest="mrna_degrad",
                  help="Parametr of mRNA degradation")

    parser.add_option("-l", "--translation", default = 0.005172, dest="translation",
                  help="Parametr of translation")

    parser.add_option("-p", "--tau_path", default = "/Users/JDima/Documents/seminars/StochKit2.0.11/tau_leaping", dest="tau_path",
                  help="**REQUIRED Path to tau leaping file")

    parser.add_option("-m", "--model", default = "/Users/JDima/PycharmProjects/gene_expression/src/model.xml", dest="model",
                  help="**REQUIRED Model file name")

    parser.add_option("-t", "--time", default = "1000", dest="time",
                  help="**REQUIRED Simulation time (i.e. run each realization from t=0 to t=time)")

    parser.add_option("-r", "--realizations", default = "1", dest="realizations",
                  help="**REQUIRED Number of realizations")

    parser.add_option("-i", "--intervals", default = "1000", dest="intervals",
                  help="Number of intervals.")
    return parser.parse_args()
