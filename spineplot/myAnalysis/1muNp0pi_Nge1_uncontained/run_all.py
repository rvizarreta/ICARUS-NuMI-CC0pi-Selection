import os
import sys

# Add spineplot to path
workspace_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(workspace_root, "../../.."))

# Import variable-specific modules
sys.path.insert(0, os.path.join(workspace_root, "dpT"))
import dpT_plots

def main(close_figs=True):
    print("=== Running dpT plots ===")
    dpT_plots.run(close_figs=close_figs)

if __name__ == '__main__':
    main()