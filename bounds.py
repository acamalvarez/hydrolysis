class Bounds:
    '''
    Lower and upper bounds for the models
    '''
    # just a and b
    ###### k2, k3, KM
    lb_ab = [0, 0]
    ub_ab = [1e4, 1e4]

    # with no inhibition
    ###### k2, k3, KM
    lb_no = [0, 0, 43.91]
    ub_no = [1e4, 1e4, 43.93]

    # substrate inhibition
    ## k2, k3, Km, Ks
    lb_s = [0, 0, 43.91, 0]
    ub_s = [1e4, 1e4, 43.93, 1e4]


    # product inhibition
    ####### k2,  k3,  Km,      Kp,         p
    lb_p = [0,   0,   117.449, 22.3939,    0]
    ub_p = [1e4, 1e4, 117.451, 22.3941,    1e4]


    # substrate and product inhibition
    ## k2, k3, Km, Kp, p, Ks
    lb_sp = [0, 0, 43.91, 0, 0, 0]
    ub_sp = [1e4, 1e4, 43.93, 1e4, 1e4, 1e4]


    # product inhibition Camilo
    ####### k2
    lb_Cam = [0, 0]
    ub_Cam = [1e4, 1e4]

    # product inhibition Tatiana
    ####### kd, k2
    lb_Tati = [0, 0]
    ub_Tati = [1e4, 1e4]

    # product inhibition Qi_p
    ####### kd, k2
    lb_Qi_p = [0, 0, 0]
    ub_Qi_p = [1e4, 1e4, 1e4]