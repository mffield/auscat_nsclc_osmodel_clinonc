function [stage_string,stage] = overallstage(ts,ns,ms,sgroup)

% This function converts a set of TNM stage features to an overall stage based on lung tumour.
stage_string='missing';
stage = 0;
if strcmp(ts,'') || any(isnan(ts)); ts='?'; end
if strcmp(ns,'')|| any(isnan(ns)); ns='?'; end
if strcmp(ms,'')|| any(isnan(ms)); ms='?'; end
if strcmp(sgroup,'')|| any(isnan(sgroup)); sgroup='?'; end

if (contains(ts,'?')||contains(ts,'X')) && (contains(ns,'?')||contains(ns,'X')) && (contains(ms,'?')||contains(ms,'X'))
    if (contains(sgroup,'?')||contains(sgroup,'X'))
        return
    else
        stage_string = sgroup;
    end
    
elseif (contains(ms,'?')||contains(ms,'X'))
    if (contains(sgroup,'?')||contains(sgroup,'X'))
        
        if contains(ns,'0')
            
            if contains(ts,'is')
                stage_string = 'Stage 0';
            elseif contains(ts,'1')
                stage_string = 'Stage IA';
            elseif contains(ts,'2a')
                stage_string = 'Stage IB';
            elseif contains(ts,'1') || contains(ts,'2a') || contains(ts,'2 S') || contains(ts,'is')
                stage_string = 'Stage I';
            elseif contains(ts,'2b')
                stage_string = 'Stage IIA';
            elseif contains(ts,'3')
                stage_string = 'Stage IIB';
            elseif contains(ts,'4')
                stage_string = 'Stage IIIA';
            else
                return
            end
        elseif contains(ns,'1')
            
            if contains(ts,'0') || contains(ts,'1') || contains(ts,'2a') || contains(ts,'2 S')
                stage_string = 'Stage IIA';
            elseif contains(ts,'2b')
                stage_string = 'Stage IIB';
            elseif contains(ts,'3') || contains(ts,'4')
                stage_string = 'Stage IIIA';
            end
            
        elseif contains(ns,'2')
            if contains(ts,'0') || contains(ts,'1') || contains(ts,'2') || contains(ts,'3') 
                stage_string = 'Stage IIIA';
            elseif contains(ts,'4')
                stage_string = 'Stage IIIB';
            end
            
        elseif contains(ns,'3')
            stage_string = 'Stage IIIB';
            
        end
    else
        stage_string = sgroup;
    end
    
    
elseif contains(ms,'1')
    stage_string = 'Stage IV';

    
elseif contains(ns,'0')

    if contains(ts,'is')
            stage_string = 'Stage 0';
    elseif contains(ts,'1')
            stage_string = 'Stage IA';
    elseif contains(ts,'2a')
            stage_string = 'Stage IB';
    elseif contains(ts,'1') || contains(ts,'2a') || contains(ts,'2 S') || contains(ts,'X') || contains(ts,'is') || contains(ts,'?')
            stage_string = 'Stage I';
    elseif contains(ts,'2b')
            stage_string = 'Stage IIA';
    elseif contains(ts,'3')
            stage_string = 'Stage IIB';
    elseif contains(ts,'4')
            stage_string = 'Stage IIIA';
    end
    
elseif contains(ns,'1')
    
    if contains(ts,'0') || contains(ts,'1') || contains(ts,'2a') || contains(ts,'2 S') || contains(ts,'X') || contains(ts,'?')
            stage_string = 'Stage IIA';
    elseif contains(ts,'2b')
            stage_string = 'Stage IIB';
    elseif contains(ts,'3') || contains(ts,'4')
            stage_string = 'Stage IIIA';     
    end
    
elseif contains(ns,'2')
    if contains(ts,'0') || contains(ts,'1') || contains(ts,'2') || contains(ts,'3') || contains(ts,'X') || contains(ts,'?')
            stage_string = 'Stage IIIA';
    elseif contains(ts,'4')
            stage_string = 'Stage IIIB';   
    end
   
elseif contains(ns,'3')
        stage_string = 'Stage IIIB';
        
elseif contains(ms,'0') && (contains(ns,'X') || contains(ns,'?') || contains(ts,'X') || contains(ts,'?'))
        if contains(sgroup,'?')
            stage_string = 'Stage Not IV';
        else
            stage_string = sgroup;
        end
elseif contains(ns,'X') || contains(ns,'?')
        stage_string = 'Stage missing';
    
elseif (contains(ms,'?') || contains(ms,'X'))
    stage_string='missing'; stage=0;
    
end

end

