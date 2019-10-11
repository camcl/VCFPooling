import os
import sys
sys.path.append(os.path.dirname(__file__))
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
print(sys.path)

# __all__ = ["alleles",
#            "pool",
#            "beagle_impute_v0.x.py",
#            "beagle_tools",
#            "bcftools",
#            "pools_plots"]
