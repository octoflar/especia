//! @file dataio.h
//! Data input and output procedures.
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#ifndef ESPECIA_DATAIO_H
#define ESPECIA_DATAIO_H

#include <iostream>
#include <valarray>

#include "base.h"

namespace especia {

    /// Reads spectroscopic data from an input stream (two-column format).
    ///
    /// @param[in,out] is The input stream.
    /// @param[out] x The wavelength data (previous content is overwritten).
    /// @param[out] y The spectral flux or uncertainty data (previous content is overwritten).
    /// @param[in] skip The number of leading lines to skip (optional).
    ///
    /// @return the input stream.
    std::istream &
    get(std::istream &is, std::valarray<real> &x, std::valarray<real> &y, natural skip = 0);

    /// Reads spectroscopic data from an input stream.
    ///
    /// @param[in,out] is The input stream.
    /// @param[out] x The wavelength data (previous content is overwritten).
    /// @param[out] y The spectral flux data (previous content is overwritten).
    /// @param[out] z The spectral flux uncertainty data (previous content is overwritten).
    /// @param[in] skip The number of leading lines to skip (optional).
    ///
    /// @return the input stream.
    std::istream &
    get(std::istream &is, std::valarray<real> &x, std::valarray<real> &y,
        std::valarray<real> &z, natural skip = 0);

    /// Writes spectroscopic data to an output stream.
    ///
    /// @param[in,out] os The output stream.
    /// @param[in] x The wavelength data.
    /// @param[in] y The spectral flux data.
    /// @param[in] z The spectral flux uncertainty data.
    ///
    /// @return the output stream.
    std::ostream &
    put(std::ostream &os, const std::valarray<real> &x, const std::valarray<real> &y,
        const std::valarray<real> &z);

}

#endif // ESPECIA_DATAIO_H
