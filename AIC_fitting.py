



def AIC_fitting(Func, x_1,x_2,y_1,y_2, x_data, y_data, y_error):

    AIC = []
    for i in range(x_1,x_2):
        for j in range(y_1,y_2):
            param, povc = curve_fit(Func, x_data[i,j], y_data[i,j], y_error[i,j]) 

            res1 = 0
            y1 = y1[i:j]
            for k in range(len(y1)):
                RSS += ((Func(x_data[i:j], *param)[k]-y1[k])**2)
            AIC.append([i,j,np.exp(4-2*np.log(RSS/len(y1))), *param, *povc])
    
    AIC = np.array(AIC)
    AICmin = np.min(AIC[:,2])

    summass = 0 
    sumerr  = 0
    probsum = 0
    for i in range(len(AIC[:,0])):
        probsum += np.exp((AICmin - AIC[i,3])/2)
        if len(AIC[0,:]) == 7: #one exponential fit
            summass += np.exp((AICmin - AIC[i,3])/2)*AIC[i,4]
            sumerr  += np.exp((AICmin - AIC[i,3])/2)*AIC[i,6]
        if len(AIC[0,:]) == 11: #two expontential fit
            summass += np.exp((AICmin - AIC[i,3])/2)*np.min([AIC[i,4], AIC[i,6]])
            sumerr  += np.exp((AICmin - AIC[i,3])/2)*np.min([AIC[i,8], AIC[i,10]])
        if len(AIC[0,:]) == 9: #exponential + c
            summass += np.exp((AICmin - AIC[i,3])/2)*AIC[i,5]
            sumerr  += np.exp((AICmin - AIC[i,3])/2)*AIC[i,8]

    mass = summass/probsum
    mass_err = sumerr/probsum
    
    fit_params = []
    for i in range(len(AIC[:,0])):
        if AIC[i,3] == AICmin:
            fit_params.append(AIC[i,:])


    return mass, mass_err, fit_params

