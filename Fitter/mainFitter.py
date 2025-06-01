from Fitter.DosFitter import *
import argparse
import yaml
import logging

def loadConfig(configPath):
  with open(configPath, 'r') as f:
    return yaml.safe_load(f)

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", type=str, default="Fitter/fitConfig.yaml", help="Path to config file")
  args = parser.parse_args()

  print("Args parsed", flush=True)

  config = loadConfig(args.config)
  print("Config read", config, flush=True)


  dosFit = DosFitter(runsDir=config["runsDir"],
                     dosExpPath=config["dosExpPath"],
                     eMax=config["eMax"])
  dosFit.fit(paramBounds=config["bounds"])

if __name__ == "__main__":
  main()
