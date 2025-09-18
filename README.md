Compile scram b -j 4

Test cmsRun Anl_*_HT_base_cfg*.py

Run crab jobs for trigger efficiency calculations crab submit Crab_Trigger*_Run2023*_v1_*cfg.py

Run condor jobs to create plots using batch_sub_trig_trk.py

Create plots using  plot_Trigger_nTrk.py
