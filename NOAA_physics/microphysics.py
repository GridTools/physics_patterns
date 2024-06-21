QCMIN = 1.e-6
TICE0 = 273.16


def moist_heat_capacity(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    c1_ice,
    c1_liq,
    c1_vap,
):
    q_liq = qliquid + qrain
    q_solid = qice + qsnow + qgraupel
    return 1.0 + qvapor * c1_vap + q_liq * c1_liq + q_solid * c1_ice


def write_with_vertical_offsets(
    qvapor,
    qliquid,
    qrain,
    qice,
    qsnow,
    qgraupel,
    cvm,
    temperature,
    delp,
    z_edge,
    z_terminal,
    z_surface,
    timestep,
    v_terminal,
    r1,
    tau_mlt,
    icpk,
    li00,
    c1_vap,
    c1_liq,
    c1_ice,
    ks,
    ke,
    is_,
    ie,
    js,
    je,
):
    '''
    This is an example of a double k-loop with assignments to multiple k-levels
    in the loop. We loop from the bottom of the atmosphere up in the k-loop and
    if there is enough ice at k to melt and precipitate we then loop back down
    in the atmosphere (the m-loop), precipitating down from k to m, and update
    the tracer ratios, temperatures, and heat capacities at both k and m.
    '''
    for i in range(is_, ie + 1):
        for j in range(js, je + 1):
            for k in range(ke - 1, ks - 1, -1):
                if v_terminal[i, j, k] < 1.0e-10:
                    continue
                if qice[i, j, k] > QCMIN:
                    for m in range(k + 1, ke + 1):
                        if z_terminal[i, j, k + 1] >= z_edge[i, j, m]:
                            break
                        if (z_terminal[i, j, k] < z_edge[i, j, m + 1]) and (
                            temperature[i, j, m] > TICE0
                        ):
                            cvm[i, j, k] = moist_heat_capacity(
                                qvapor[i, j, k],
                                qliquid[i, j, k],
                                qrain[i, j, k],
                                qice[i, j, k],
                                qsnow[i, j, k],
                                qgraupel[i, j, k],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            cvm[i, j, m] = moist_heat_capacity(
                                qvapor[i, j, m],
                                qliquid[i, j, m],
                                qrain[i, j, m],
                                qice[i, j, m],
                                qsnow[i, j, m],
                                qgraupel[i, j, m],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            dtime = min(
                                timestep,
                                (z_edge[i, j, m] - z_edge[i, j, m + 1])
                                / v_terminal[i, j, k],
                            )
                            dtime = min(1.0, dtime / tau_mlt)
                            sink = min(
                                qice[i, j, k] * delp[i, j, k] / delp[i, j, m],
                                dtime
                                * (temperature[i, j, m] - TICE0)
                                / icpk[i, j, m],
                            )
                            qice[i, j, k] -= sink * (
                                delp[i, j, m] / delp[i, j, k]
                            )
                            if z_terminal[i, j, k] < z_surface[i, j]:
                                r1[i, j] += sink * delp[i, j, m]
                            else:
                                qrain[i, j, m] += sink

                            cvm_tmp = moist_heat_capacity(
                                qvapor[i, j, k],
                                qliquid[i, j, k],
                                qrain[i, j, k],
                                qice[i, j, k],
                                qsnow[i, j, k],
                                qgraupel[i, j, k],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            temperature[i, j, k] = (
                                temperature[i, j, k] * cvm[i, j, k]
                                - li00 * sink * delp[i, j, m] / delp[i, j, k]
                            ) / cvm_tmp
                            cvm_tmp = moist_heat_capacity(
                                qvapor[i, j, m],
                                qliquid[i, j, m],
                                qrain[i, j, m],
                                qice[i, j, m],
                                qsnow[i, j, m],
                                qgraupel[i, j, m],
                                c1_ice,
                                c1_liq,
                                c1_vap,
                            )
                            temperature[i, j, m] = (
                                temperature[i, j, m] * cvm[i, j, m]
                            ) / cvm_tmp
                        if qice[i, j, k] < QCMIN:
                            break
