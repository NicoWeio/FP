import contextlib
import numpy as np
import os
import pint
import scipy as sp
import scipy.stats
from uncertainties import unumpy, ufloat, UFloat

### Umgebungsvariablen ###
BUILD = bool(os.getenv('BUILD'))
print("build:", BUILD)
# TODO: BUILD steuert gerade verschiedene Dinge:
# 1. sollen Plots gespeichert werden?
# 2. sollen Plots mit LaTeX generiert werden?

if not BUILD:
    # nutze Standard-Backend, nicht LaTeX
    os.environ['MATPLOTLIBRC'] = ''

### Funktionen ###


def linregress(x, y):
    r = sp.stats.linregress(x.m, y.m)
    return (
        ufloat(r.slope, r.stderr) * (y.units / x.units),
        ufloat(r.intercept, r.intercept_stderr) * y.units
    )


def curve_fit(fit_fn, x, y, **kwargs):
    """Wrapper for `sp.optimize.curve_fit` that returns `ufloat`s"""
    if kwargs.get('bounds', False) is None:
        del kwargs['bounds']

    params, pcov = sp.optimize.curve_fit(fit_fn, x, y, **kwargs)
    param_errors = np.sqrt(np.diag(pcov))
    return tuple(ufloat(p, e) for p, e in zip(params, param_errors))


def pint_curve_fit(fit_fn, x, y, param_units, bounds=None, p0=None):
    """Wrapper for `sp.optimize.curve_fit` that accepts pint units."""
    # TODO: Abweichung mittels sigma-Parameter berücksichtigen

    def convert(value, unit):
        # if value is None:
        #     return None
        return value.to(unit).m

    # if x.units == ureg.deg or y.units == ureg.deg:
    #     print("⚠️ You are trying to pass degrees to a function that probably expects radians!")
    #     # Since fit_fn only receives magnitudes, degrees won't behave as expected in e.g. np.sin(…).
    #     # TODO: Best solution would be to write a wrapper for fit_fn

    if bounds is not None:
        # first, [(lower, upper)]
        bounds = [
            (b[0].to(u).m, b[1].to(u).m) if (b is not None) else (-np.infty, np.infty)
            for b, u in zip(bounds, param_units)
        ]
        # then, ([lower], [uppper])
        bounds = tuple(zip(*bounds))

    if p0 is not None:
        p0 = [convert(p, unit) for p, unit in zip(p0, param_units, strict=True)]

    u_params = curve_fit(fit_fn, x.m, y.m, p0=p0, bounds=bounds)
    pint_params = tuple(p * u for p, u in zip(u_params, param_units))

    try:
        pint_params_nominal = tuple(p.n * u for p, u in zip(u_params, param_units))
        test_val = fit_fn(x, *pint_params_nominal)
    except:
        raise Exception("Could not test fit_fn")
    if test_val.units != y.units:
        raise Exception(
            f"Wrong param_units – fit_fn(x[0], *fit_params_nominal) returns '{test_val.units}' instead of '{y.units}'")
    return pint_params


def pint_polyfit(x, y, deg):
    params, covariance_matrix = np.polyfit(x.m, y.m, deg=deg, cov=True)
    errors = np.sqrt(np.diag(covariance_matrix))
    return [ufloat(param, error) * y.units / x.units**(deg-i) for i, (param, error) in enumerate(zip(params, errors))]


def pintify(list):
    assert len(list) > 0
    units = list[0].units
    # assert all(e.units == units for e in list)
    return [e.to(units).m for e in list] * units


def uarray(nominal_values, std_devs):
    # assert len(nominal_values) == len(std_devs)
    units = nominal_values.units
    return unumpy.uarray(nominal_values.to(units).m, std_devs.to(units).m) * units


def nominal_values(list):
    assert isinstance(list, pint.Quantity)
    units = list.units
    return [e.m.n for e in list] * units


def std_devs(list):
    assert isinstance(list, pint.Quantity)
    units = list.units
    return [e.m.s for e in list] * units


def nominal_value(v):
    units = v.units
    return v.m.n * units


def linspace(start, end, num=50):
    return np.linspace(start.m, end.to(start.units).m, num=num) * start.units


def ufloat_from_list(vals):
    if isinstance(vals, pint.Quantity):
        return ufloat(np.mean(vals.m), sp.stats.sem(vals.m)) * vals.units
    else:
        return ufloat(np.mean(vals), sp.stats.sem(vals))


def fmt_abs_err(o, r, precise=False, show_uncertainty=True):
    o_abs_err = o - r
    if isinstance(o_abs_err, UFloat) and not show_uncertainty:
        o_abs_err = o_abs_err.n
    return f'{o_abs_err:.2f}'


def fmt_rel_err_percent(o, r, precise=False, show_uncertainty=True):
    o_rel_err = ((o - r)/r).to('dimensionless').m
    if isinstance(o_rel_err, UFloat) and not show_uncertainty:
        o_rel_err = o_rel_err.n
    return f'{o_rel_err:.2%}'


def fmt_err(o, r, precise=False):
    return f'{fmt_rel_err_percent(o, r, precise)} | {fmt_abs_err(o, r, precise)}'


def fmt_compare_to_ref(o, r, name=None, precision=2, unit=None):
    my_o = o.to(unit) if unit else o
    my_r = r.to(unit) if unit else r
    my_name = f'{name}:\n' if name else ''
    return my_name + (
        f'- ist: {my_o:.{precision}f}\n'
        f'- soll: {my_r:.{precision}f}\n'
        f'- abs. Abweichung: {fmt_abs_err(my_o, my_r)}\n'
        f'- rel. Abweichung: {fmt_rel_err_percent(o, r)}'
    )


def pint_concat(*lists):
    units = lists[0].units
    out = []
    for l in lists:
        vals = l.to(units).m
        out.extend(vals)
    out *= units
    return out
    # return [*[l.to(units).m for l in lists]] * units


# Entfernt Wertepaare(/-tupel…), die NaNs enthalten
def remove_nans(*inputs):
    # return tuple(zip(*[input_tuple for input_tuple in zip(*inputs) if not any(np.isnan(v) for v in input_tuple)]))
    return (pintify(x) for x in zip(*[input_tuple for input_tuple in zip(*inputs) if not any(np.isnan(v) for v in input_tuple)]))


# Hilfreich, um z.B. eine Gerade zu plotten…
def bounds(vals):
    return pintify([min(vals), max(vals)])


def errorbar(plt, x, y, **kwargs):
    def get_n_s(vals):
        """trennt vals in nominal_values und std_devs und gibt sie ohne Einheit zurück"""
        try:
            return nominal_values(vals).m, std_devs(vals).m
        except AttributeError:  # scheinbar keine Unsicherheiten angegeben
            return vals.m, None

    x_n, x_s = get_n_s(x)
    y_n, y_s = get_n_s(y)

    return plt.errorbar(x_n, y_n, xerr=x_s, yerr=y_s, **kwargs)


@contextlib.contextmanager
def plot_context(plt, xunits, yunits, xname=None, yname=None):
    """Context-Manager zum einfacheren Plotten von einheitenbehafteten Daten."""

    # Parse Einheiten. Das funktioniert auch, wenn sie bereits Instanzen von pint.Unit sind.
    xunits = pint.Unit(xunits)
    yunits = pint.Unit(yunits)

    class MyPlotter:
        """Hilfsklasse, die einen Teil von `matplotlib.pyplot` nachahmt und ergänzt."""

        def plot(self, *args, show_xerr=True, show_yerr=True, **kwargs):
            x, y, *more_args = args

            # Bringe x und y in die vorgegebene Einheit und entferne sie anschließend.
            if isinstance(x, pint.Quantity):
                x = x.to(xunits).m
            if isinstance(y, pint.Quantity):
                y = y.to(yunits).m

            # Trenne x und y in nominal_values und std_devs auf.
            def get_n_s(vals):
                n, s = unumpy.nominal_values(vals), unumpy.std_devs(vals)
                # Falls alle Unsicherheiten 0 sind, wird s auf None gesetzt.
                if not s.any():
                    s = None
                return n, s

            x_n, x_s = get_n_s(x)
            y_n, y_s = get_n_s(y)

            # xerr und yerr können explizit deaktiviert werden.
            if not show_xerr:
                x_s = None
            if not show_yerr:
                y_s = None

            # TODO Entwurf: Stelle Unsicherheit mit Farbfüllung statt Fehlerbalken dar.
            if show_yerr == 'fill':
                plt.fill_between(x_n, y_n - y_s, y_n + y_s, alpha=0.2)
                y_s = None

            # Plotte die Daten mit `errorbar`, falls Unsicherheiten angegeben wurden, sonst mit `plot`.
            if (x_s is not None) or (y_s is not None):
                # return plt.errorbar(x_n, y_n, *more_args, xerr=x_s, yerr=y_s, **kwargs)
                assert not more_args, 'more_args → Fehler, versuche stattdessen kwargs 🤷🏼‍♂️'
                return plt.errorbar(x_n, y_n, xerr=x_s, yerr=y_s, **kwargs)
            else:
                return plt.plot(x_n, y_n, *more_args, **kwargs)

    # automatische Labels
    def fmt_label(name, units):
        if units == pint.Unit('dimensionless'):
            return f"${name}$"
        else:
            if BUILD:
                return f"${name}" + r" \mathbin{/} " + f"{units:Lx}$"
            else:
                return f"${name}$" + r" / " + f"{units}"

    if xname:
        plt.xlabel(fmt_label(xname, xunits))
    if yname:
        plt.ylabel(fmt_label(yname, yunits))

    yield MyPlotter()
