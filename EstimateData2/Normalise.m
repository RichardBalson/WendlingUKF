function [Normalised_Data, Multiplier] = Normalise(Data,Background_amp,fs,Period)

Samples = length(Data);
Duration = Samples/fs;
Period_samples = fs*Period;
Windows = floor(Duration/Period);
MaxD = zeros(Windows,1);
for j =1:Windows    
    MaxD(j) = max(abs(Data((j-1)*Period_samples+1:j*Period_samples)));
end
Multiplier = Background_amp/min(MaxD);
Normalised_Data = Data* Multiplier;

